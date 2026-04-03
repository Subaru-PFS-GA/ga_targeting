import os
from glob import glob
import commentjson as json
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u

from pfs.datamodel import TargetType, PfsDesign
from pfs.ga.common.scripts import Progress
import pfs.ga.common.util.astro as astro

from ...config.netflow import NetflowConfig
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_SECTORS
from ..targetingscript import TargetingScript

from ...setup_logger import logger

class CatalogScript(TargetingScript, Progress):
    """
    Command-line script to generate a catalog with the available photometry
    that includes all targets.
    """

    def __init__(self):
        super().__init__()

    def _add_args(self):
        TargetingScript._add_args(self)
        Progress._add_args(self)

        self.add_arg('--dsph', type=str, choices=DSPH_FIELDS.keys(), help='Name of a predefined dSph target.')
        self.add_arg('--m31', type=str, choices=M31_SECTORS, help='Name of a predefined M31 sector.')

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')
        self.add_arg('--in', type=str, required=True, help='Input directory containing the netflow output.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')    

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self._indir = self.get_arg('in', args, self._indir)
        self._outdir = self.get_arg('out', args, self._outdir)

    def _init_from_args(self, args):
        TargetingScript._init_from_args(self, args)
        Progress._init_from_args(self, args)

        self._create_field_from_args(args)
        self._create_config_from_field()

        # Load the configuration template files and merge with the default config
        self._load_config_files(args)

    def prepare(self):
        """
        Prepare the script for execution.
        """
        super().prepare()

        # Create the output directory
        if os.path.isdir(self._outdir):
            raise Exception(f'Output directory already exists: {self._outdir}')
        else:
            logger.debug(f'Creating output directory `{self._outdir}`')
            os.makedirs(self._outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self._outdir, os.path.basename(self.log_file))

    def run(self):
        # Load the netflow targets, these are the primary coordinates
        targets = self.load_netflow_targets()

        # Exclude sky targets
        targets = targets[targets['prefix'] != 'sky']

        # DEBUG: Get a random sample of a 1000 targets for testing
        # targets = targets.sample(1000, random_state=42)

        # Load the original target list files
        target_lists = self.load_preprocessed_target_lists()
        
        # Reindex targets lists by __target_idx for easier lookup; exclude sky targets
        target_lists = {
            k: target_list.data.set_index('__target_idx', verify_integrity=True)
            for k, target_list in target_lists.items()
            if k != 'sky'
        }

        # Load the final, merged assignments list
        assignments = self._load_fiber_assignments_all_list()
        assignments.set_index('__target_idx', inplace=True)

        # Count the number of rows in assignments for each target index then
        # join it with targets on __target_idx to get the number of assignments for each target
        assigned = assignments.groupby('__target_idx').size().rename('assigned')
        targets = targets.join(assigned, how='left')
        targets['assigned'] = targets['assigned'].fillna(0).astype(int)

        logger.info(f'Joined targets with assignments, starting to build catalog with {targets.shape[0]} targets.')

        results = []
        id_prefix = self._get_id_prefix()

        dtypes = {
            'target_idx': pd.Int64Dtype(),
            'objid': pd.Int64Dtype(),
            'primary_catalog': pd.StringDtype(),
            'primary_targetid': pd.Int64Dtype(),
            'target_type': pd.StringDtype(),
            'RA': float,
            'Dec': float,
            'pmra': float,
            'pmdec': float,
            'parallax': float,
            'rv': float,
            'epoch': float,
            'penalty': float,
            'exp_time': float,
            'priority': pd.Int16Dtype(),
            'target_class': pd.StringDtype(),
            'done_visits': pd.Int16Dtype(),
            'req_visits': pd.Int16Dtype(),
            'assigned': pd.Int16Dtype()
        }

        # For each row in targets, find the corresponding items in the targets lists
        q = 0
        for target_idx, row in tqdm(targets.iterrows(), total=targets.shape[0]):

            q += 1
            if self.top == q:
                logger.info(f'Reached top={self.top} limit, stopping after {q} targets.')
                break

            if row['prefix'] == 'sci':
                objid = id_prefix | target_idx
            else:
                objid = row['targetid']

            res = dict(
                target_idx = target_idx,
                objid = objid,
                primary_catalog = row['key'],
                primary_targetid = row['targetid'],
                target_type = row['prefix'],
                RA = row['RA'],
                Dec = row['Dec'],
                pmra = row['pmra'],
                pmdec = row['pmdec'],
                parallax = row['parallax'],
                rv = row['rv'],
                epoch = row['epoch'],
                penalty = row['penalty'],
                exp_time = row['exp_time'],
                priority = row['priority'],
                target_class = row['class'],
                done_visits = row['done_visits'],
                req_visits = row['req_visits'],
                assigned = row['assigned']
            )

            results.append(res)

            c1 = SkyCoord(row['RA']*u.deg, row['Dec']*u.deg)

            # Find the target in any of the target lists
            for key, target_list in target_lists.items():
                if target_idx not in target_list.index:
                    continue
                
                other = target_list.loc[target_idx]

                c2 = SkyCoord(other['RA']*u.deg, other['Dec']*u.deg)

                # Verify that the coordinates match
                sep = c1.separation(c2).arcsecond
                if sep > 0.1:
                    continue

                # This is considered a good match so we process this item further
                res[f'targetid_{key}'] = other['targetid']
                res[f'RA_{key}'] = other['RA']
                res[f'Dec_{key}'] = other['Dec']

                dtypes[f'targetid_{key}'] = pd.Int64Dtype()
                dtypes[f'RA_{key}'] = float
                dtypes[f'Dec_{key}'] = float
    
                # Get the available photometry for this target from the target list
                photometry = self._config.targets[key].photometry
                
                if photometry.filters is not None:
                    self.__append_magnitudes_filters(res, other, photometry)

                if photometry.bands is not None:
                    self.__append_magnitudes_bands(res, other, photometry)

        # Get the union of all columns in results
        all_columns = []
        for res in results:
            for col in res.keys():
                if col not in all_columns:
                    all_columns.append(col)

        # Combine results into a dict of lists, filling missing values with NaN
        combined_results = {col: [] for col in all_columns}
        for res in results:
            for col in all_columns:
                combined_results[col].append(res.get(col, np.nan))

        # Convert each column in combined_results to a pd.Series of the appropriate data type
        combined_results = {
            col: pd.Series(combined_results[col], dtype=dtypes.get(col, object))
            for col in all_columns
        }

        df = pd.DataFrame(combined_results)
        
        filename = os.path.join(self._outdir, 'catalog.feather')
        df.to_feather(filename)
        logger.info(f'Wrote catalog with {df.shape[0]} targets to {filename}')

    def __append_magnitudes_filters(self, res, other, photometry):
        # Column names for each filter
        for filter, filter_config in photometry.filters.items():
            for kind in ['', 'psf_', 'fiber_', 'total_']:
                mag_col = filter_config.get(f'{kind}mag')
                mag_err_col = filter_config.get(f'{kind}mag_err')
                flux_col = filter_config.get(f'{kind}flux')
                flux_err_col = filter_config.get(f'{kind}flux_err')

                if self.__append_magnitude(res, other, filter, mag_col, mag_err_col, flux_col, flux_err_col):
                    # Magnitude found for this kind, don't check the rest
                    break

    def __append_magnitudes_bands(self, res, other, photometry):
        # Column names for each band, look up filter name in column
        for band, band_config in photometry.bands.items():
            filter = other[band_config['filter']]        
            for kind in ['', 'psf_', 'fiber_', 'total_']:
                mag_col = band_config.get(f'{kind}mag')
                mag_err_col = band_config.get(f'{kind}mag_err')
                flux_col = band_config.get(f'{kind}flux')
                flux_err_col = band_config.get(f'{kind}flux_err')

                if self.__append_magnitude(res, other, filter, mag_col, mag_err_col, flux_col, flux_err_col):
                    # Magnitude found for this kind, don't check the rest
                    break

    def __append_magnitude(self, res, other, filter, mag_col, mag_err_col, flux_col, flux_err_col):
        cname_mag = f'{filter}_mag'
        cname_mag_err = f'{filter}_mag_err'

        # If the magnitude column is specified, we just use it
        if mag_col is not None:
            if cname_mag not in res:
                res[cname_mag] = other[mag_col]

                if mag_err_col is not None:
                    res[cname_mag_err] = other[mag_err_col]

            # Magnitude found for filter
            return True

        # If the flux column is specified, we can convert it to magnitude
        if flux_col is not None:
            if cname_mag not in res:
                flux = other[flux_col]
                if flux_err_col is not None:
                    flux_err = other[flux_err_col]
                else:
                    flux_err = np.nan

                mag, mag_err = astro.nJy_to_ABmag(flux, flux_err=flux_err)

                res[cname_mag] = mag
                res[cname_mag_err] = mag_err

            # Magnitude found for filter
            return True
        
        # No magnitude or flux column found for this filter
        return False

    def load_netflow_targets(self):
        targets = None

        for dir in self._indir if isinstance(self._indir, list) else [self._indir]:
            fn = os.path.join(dir, 'netflow_targets.feather')
            df = pd.read_feather(fn)

            if targets is None:
                targets = df
            else:
                targets = pd.concat([targets, df])

        return targets

    # TODO: merge function with the one in netflowscript.py
    def load_preprocessed_target_lists(self):
        target_lists = {}
        for k, target_list_config in self._config.targets.items():
            # Load the preprocessed target lists from the input directory
            # These will be later update with the number of scheduled visits
            # and saved to the output directory
            path = self._get_preprocessed_target_list_path(k, self._indir)
            target_lists[k] = self._load_preprocessed_target_list(k, path)

        return target_lists

    # TODO: merge function with the one in exportscript.py
    def _load_fiber_assignments_all_list(self):
        # Load the list of fiber assigments from each input directory
        fiber_assignemnts_all_list = None
        for dir in self._indir if isinstance(self._indir, list) else [self._indir]:
            fn = self._get_fiber_assignments_all_path(dir)
            df = pd.read_feather(fn)

            if fiber_assignemnts_all_list is None:
                fiber_assignemnts_all_list = df
            else:
                fiber_assignemnts_all_list = pd.concat([fiber_assignemnts_all_list, df])

        return fiber_assignemnts_all_list

if __name__ == '__main__':
    script = CatalogScript()
    script.execute()