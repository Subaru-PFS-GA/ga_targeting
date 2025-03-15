import os
import commentjson as json
from datetime import datetime, timedelta, tzinfo
import numpy as np
import pandas as pd
from pandas import Float32Dtype, Float64Dtype, Int32Dtype, Int64Dtype
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

from pfs.datamodel import TargetType, FiberStatus, PfsDesign

import pfs.ga.targeting
from ...config.netflow import NetflowConfig
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_FIELDS
from ...instrument import SubaruPFI
from ...io import ObservationSerializer
from ...netflow import Netflow, Design
from ...util.args import *
from ...util.astro import *
from ...util.pandas import *
from ..targetingscript import TargetingScript

from ...setup_logger import logger

class ImportScript(TargetingScript):
    """
    Command-line script to import target lists into a common format, ready for
    running netflow in stages.

    This scripts import the target lists, converts them into a uniform format,
    performs the necessary coordinate transformations and
    """

    # TODO: see what can be merged from ImportScript, ExportScript and NetflowScript into TargetingScript

    def __init__(self):
        super().__init__()

        self.__exp_time = None
        self.__obs_time = None
        self.__xmatch_rad = 2           # arcseconds
        self.__skip_notebooks = False

    def _add_args(self):
        super()._add_args()

        self.add_arg('--dsph', type=str, choices=DSPH_FIELDS.keys(), help='Name of a predefined dSph target.')
        self.add_arg('--m31', type=str, choices=M31_FIELDS, help='Name of a predefined M31 field.')

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')
        self.add_arg('--exp-time', type=float, help='Exposure time per visit, in seconds.')
        self.add_arg('--obs-time', type=str, help='Observation time in ISO format in UTC.')        
        self.add_arg('--xmatch-rad', type=float, help='Cross-match radius in arcseconds.')
        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')    

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self._outdir = self.get_arg('out', args, self._outdir)

    def _init_from_args(self, args):
        super()._init_from_args(args)

        # If a field is specified, load its default configuration      
        self._create_field_from_args(args)
        self._create_config_from_field()

        # Load the configuration template files and merge with the default config
        self._load_config_files(args)

        self.__exp_time = self.get_arg('exp_time', args, self.__exp_time)
        self.__obs_time = self.get_arg('obs_time', args, self.__obs_time)
        self.__xmatch_rad = self.get_arg('xmatch_rad', args, self.__xmatch_rad)
        self.__skip_notebooks = self.get_arg('skip_notebooks', args, self.__skip_notebooks)

        # Override the configuration with the command-line arguments
        if self.__exp_time is not None:
            self._config.field.exp_time = self.__exp_time

        if self.__obs_time is not None:
            self.__obs_time = datetime.fromisoformat(self.__obs_time)
            self._config.field.obs_time = self.__obs_time

    def prepare(self):
        super().prepare()

        # Create the output directory
        if os.path.isdir(self._outdir):
            raise Exception(f'Output directory already exists: {self._outdir}')
        else:
            os.makedirs(self._outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self._outdir, os.path.basename(self.log_file))

    def run(self):

        # Load instrument calibration data
        instrument = self._create_instrument()

        # Generate the list of pointings
        pointings = self._generate_pointings()

        # Create the netflow object
        netflow = Netflow(
            f'{self._config.field.name}',
            instrument,
            pointings,
            workdir = self._outdir,
            netflow_options = self._config.netflow_options,
            solver_options = self._config.gurobi_options,
            debug_options = self._config.debug_options)

        # Load the original target list files or, if resuming the processing,
        # load the preprocessed target lists from the output directory.
        target_lists = self.load_source_target_lists()

        # Look for required columns, duplicates, etc.
        self.__validate_source_target_lists(target_lists)

        # Transform the coordinates to a common epoch and ICRS
        self.__transform_source_target_lists(target_lists)

        # Append the source target lists to the netflow object
        # This will perform the cross-matching and generate the target indices
        # that are written back into the target list DataFrames
        # This target_idx will be unique
        self.__append_source_target_lists(netflow, target_lists)

        # Because the target lists have changed (__key and __target_idx added)
        # we save them again at this point
        self._save_preprocessed_target_lists(target_lists, self._outdir)

        # Save the netflow targets, by default to the output directory.
        netflow.validate_targets()
        netflow.save_targets()

        # Execute the evaluation notebooks
        if not self.__skip_notebooks:
            for notebook in ['targets']:
                logger.info(f'Executing evaluation notebook `{notebook}`...')
                notebook_path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), f'scripts/import/notebooks/{notebook}.ipynb')
                parameters = {
                    'DEBUG': False,
                    'CONFIG_FILE': self._get_output_config_path(),
                    'OUTPUT_PATH': self._outdir,
                }
                self._execute_notebook(notebook_path, parameters, self._outdir)
    
    def load_source_target_lists(self):
        target_lists = {}
        for k, target_list_config in self._config.targets.items():
            target_list = self.load_target_list(k, target_list_config)
            target_list.name = k
            
            self.__calculate_target_list_flux(k, target_list)
            self.__append_target_list_extra_columns(k, target_list)
            self.__validate_target_list(k, target_list_config, target_list)
            
            target_lists[k] = target_list

        return target_lists

    def __calculate_target_list_flux(self, key, target_list):
        # TODO: this should go under the Observation class, maybe

        # Calculate the missing flux columns
        # If no psf, fiber, total, etc. flux/magnitude is available, just copy the values

        # If the entry in the config is a dict we assume that filter names come from the column
        # names and the same filters are used for the entire target list. If the config entry is
        # a list, we assume that the filter names are stored in a filter column. This latter is
        # typical for flux standard lists.

        photometry = self._config.targets[key].photometry

        if photometry is not None and photometry.filters is not None:
            return self.__calculate_target_list_flux_filters(key, target_list)
        elif photometry is not None and  photometry.bands is not None:
            return self.__calculate_target_list_flux_bands(key, target_list)
        else:
            # No fluxes or magnitudes available (sky)
            pass

    def __calculate_target_list_flux_filters(self, key, target_list):
        """
        Calculate the missing flux values from other columns that are available
        """

        target_list.calculate_flux_filters(self._config.targets[key].photometry.filters)

    def __calculate_target_list_flux_bands(self, key, target_list):
        """
        Calculate the missing flux values from other columns that are available
        """

        target_list.calculate_flux_bands(self._config.targets[key].photometry.bands)

    def __append_target_list_extra_columns(self, key, target_list):
        """
        Join the assignments with the target lists to append the fluxes that will
        be included in the design files.
        """

        logger.info(f'Appending extra columns to target list `{key}`.')

        # Substitute these fields into column patterns
        field_info = dict(
            name=self._config.field.key,               # Field short name, e.g. 'umi'
            key=key,                                    # Target list key, e.g. 'dsph' or 'fluxstd'
            obs_time=self._config.field.obs_time,
            resolution=self._config.field.resolution,
            prefix=self._config.targets[key].prefix,
            catid=self._config.targets[key].catid,
        )

        def set_extra_column_constant(name, dtype, default_value, config_value):
            if config_value is not None:
                target_list.append_column(name, config_value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using the config value {config_value}.')
            elif name not in target_list.columns and default_value is not None:
                target_list.append_column(name, default_value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using the default value {default_value}.')
            elif name in target_list.columns:
                target_list.set_column_dtype(name, dtype)
                logger.debug(f'Updated dtype of column `{name}` in target list `{key}`.')

        def set_extra_column_lambda(name, dtype, default_value, func, args):
            if func is not None:
                if isinstance(args, list):
                    value = target_list.data[args].apply(func, axis=1)
                else:
                    value = target_list.data[args].apply(func)
                target_list.append_column(name, value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using a lambda function.')
            elif name not in target_list.columns and default_value is not None:
                target_list.append_column(name, default_value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using the default value {default_value}.')
            elif name in target_list.columns:
                target_list.set_column_dtype(name, dtype)
                logger.debug(f'Updated dtype of column `{name}` in target list `{key}`.')

        def set_extra_column_pattern(name, dtype, default_value, config_pattern):
            if config_pattern is not None:
                value = config_pattern.format(**field_info)
                if '{' in value:
                    value = target_list.data[['targetid']].apply(lambda row: value.format(**row), axis=1)
                target_list.append_column(name, value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using a config pattern.')
            elif name not in target_list.columns and default_value is not None:
                target_list.append_column(name, default_value, dtype)
                logger.debug(f'Added column `{name}` to target list `{key}` using the default value {default_value}.')
            elif name in target_list.columns:
                target_list.set_column_dtype(name, dtype)
                logger.debug(f'Updated dtype of column `{name}` in target list `{key}`.')

        # Verify extra columns config, only one type of pattern is allowed
        if self._config.targets[key].extra_columns is not None:
            for c, config in self._config.targets[key].extra_columns.items():
                patterns = [config.constant, config.lambda_func, config.pattern]
                if np.sum([p is not None for p in patterns]) > 1:
                    raise Exception(f'Only one type of pattern is allowed for extra column `{c}`.')
                
            # The order of column generation is: constants, lambda, pattern
            for c, config in self._config.targets[key].extra_columns.items():
                if config.constant is not None:
                    set_extra_column_constant(c, config.dtype, None, config.constant)
            
            for c, config in self._config.targets[key].extra_columns.items():
                if config.lambda_func is not None:
                    set_extra_column_lambda(c, config.dtype, None, config.lambda_func, config.lambda_args)

            for c, config in self._config.targets[key].extra_columns.items():
                if config.pattern is not None:
                    set_extra_column_pattern(c, config.dtype, None, config.pattern)

        # If constants are defined for these columns in the config, assign them
        epoch = normalize_epoch(self._config.targets[key].epoch)
        if epoch is not None:
            epoch = epoch.jyear
            target_list.epoch = epoch
            set_extra_column_constant('epoch', 'string', 'J2016.0', f'J{epoch:0.2f}')
        else:
            set_extra_column_constant('epoch', 'string', 'J2016.0', None)
           
        set_extra_column_constant('catid', pd.Int32Dtype(), -1, self._config.targets[key].catid)
        set_extra_column_constant('priority', pd.Int32Dtype(), None, self._config.targets[key].priority)
        set_extra_column_constant('exp_time', pd.Int32Dtype(), None, self._config.targets[key].exp_time)

        # Create required columns with default values
        if 'targetid' not in target_list.columns:
            target_list.append_column('targetid', target_list.data.index, pd.Int64Dtype())
            logger.warning(f'Automatically generated column `targetid` in target list `{key}`.')

        if 'proposalid' not in target_list.columns:
            target_list.append_column('proposalid', '', 'string')

        if 'obcode' not in target_list.columns:
            target_list.append_column('obcode', target_list.data['targetid'], 'string')

        # TODO: these must be calculated from the coordinates
        if 'tract' not in target_list.columns:
            target_list.append_column('tract', 0, pd.Int32Dtype())
        else:
            target_list.set_column_dtype('tract', pd.Int32Dtype())
        
        if 'patch' not in target_list.columns:
            target_list.append_column('patch', '0,0', 'string')
        else:
            target_list.set_column_dtype('patch', 'string')

        logger.info(f'Successfully added extra columns to target list `{key}`.')

    def __validate_target_list(self, key, target_list_config, target_list):
        # Verify that columns are present

        logger.info(f'Validating target list `{key}`.')

        must = ['targetid', 'RA', 'Dec']
        if target_list_config.prefix == 'sci':
            must.extend(['priority', 'exp_time'])

        for col in must:
            if col not in target_list.data.columns:
                raise Exception(f'Missing column "{col}" in target list `{key}`.')

        should = []
        if target_list_config.prefix in ['sci', 'cal']:
            should.extend(['pmra', 'pmdec', 'parallax', 'rv', 'epoch'])
            
        for col in should:
            if col not in target_list.data.columns:
                logger.warning(f'Missing column "{col}" in target list `{key}`.')

        # Verify that targetid is unique and print duplicates
        duplicates = target_list.data[target_list.data.duplicated('targetid', keep=False)]
        if not duplicates.empty:
            logger.warning(f'Duplicate targetids found in target list `{key}`.')
            logger.warning(list(duplicates['targetid']))

        # TODO: verify that all columns are present that are in the filter list
        #       epoch, catid, proposalid, etc.

        # Verify that science targets have a priority and exposure time set
        if target_list_config.prefix == 'sci':
            if np.any(target_list.data['priority'].isna()):
                raise Exception(f'Missing priority in science target list `{key}`.')
            if np.any(target_list.data['exp_time'].isna()):
                raise Exception(f'Missing exposure time in science target list `{key}`.')
            if np.any(target_list.data['exp_time'] <= 0):
                raise Exception(f'Invalid exposure time in science target list `{key}`.')

        logger.info(f'Successfully validated target list `{key}`.')

    def __validate_source_target_lists(self, target_lists):
        # TODO: make sure all targetids are unique across all target lists

        required_columns = {
            'sci': [ 'epoch', 'proposalid', 'tract', 'patch', 'catid', 'obcode' ],
            'cal': [ 'epoch', ],
            'sky': [ 'epoch', ],
            'ag': [ 'epoch', 'parallax', 'pmra', 'pmdec' ]
        }

        for key, target_list in target_lists.items():
            if self._config.targets[key].prefix not in ['sci', 'cal', 'sky', 'ag']:
                raise Exception(f'Invalid target list prefix `{self._config.targets[key].prefix}` for target list `{key}`.')

            for col in required_columns[self._config.targets[key].prefix]:
                if col not in target_list.columns:
                    raise Exception(f'Missing required column "{col}" in target list `{key}`.')

    def __transform_source_target_lists(self, target_lists):
        # Transform coordinates and apply proper motion, when necessary, to convert each catalog
        # to ICRS and the same epoch

        for key, target_list in target_lists.items():    
            self.__transform_target_coords(target_list)

    def __transform_target_coords(self, target_list):
        # Transform coordinates and apply proper motion, when necessary, to convert each catalog
        # to ICRS and the same epoch.
        # Note, that `apply_space_motion` requires `radial_velocity` to be set when working in ICRS

        # TODO: move this to the Catalog object from here?

        target_epoch = normalize_epoch(self._config.netflow_options.epoch)
        
        logger.info(f'Making sure target list `{target_list.name}` is in ICRS and at epoch J{target_epoch.jyear:0.2f}.')
        target_list.transform_coords(target_frame='icrs', target_epoch=target_epoch)

    def __append_source_target_lists(self, netflow, target_lists):
        # Add flux standards first because those are the less numerous. This way
        # we always cross-match the science targets with the flux standards and
        # can throw away stars that are already in the flux standard list.
        # Sky positions are added at the end to avoid cross-matching

        for prefix in ['cal', 'sci']:
            for k, target_list in target_lists.items():
                if self._config.targets[k].prefix == prefix:
                    # Remember the key of the target list so that it will appear in the
                    # final assigned targets's list
                    pd_append_column(target_list.data, '__key', k, pd.StringDtype())

                    # Cross-match the target list with the target lists already added to netflow
                    # idx1 is the index into the current target list
                    # idx2 is the index of the matches in the netflow target cache
                    idx1, idx2, _, mask = self.__cross_match_source_target_list(netflow, target_list)

                    # Make room for target_idx in the target list for joining in the magnitudes later
                    pd_append_column(target_list.data, '__target_idx', pd.NA, pd.Int64Dtype())
                    
                    if mask is not None:
                        # Update the target list where there is a match
                        netflow.update_targets(target_list, prefix, idx1, idx2)
                        target_list.data.loc[mask, '__target_idx'] = idx2

                        # Add targets which didn't have a match
                        idx = netflow.append_targets(target_list, prefix=prefix, mask=~mask)
                        target_list.data.loc[~mask, '__target_idx'] = idx
                    else:
                        # Add all targets as this is the first target list being processed
                        idx = netflow.append_targets(target_list, prefix=prefix)
                        target_list.data['__target_idx'] = idx
                                
        for k, target_list in target_lists.items():
            if self._config.targets[k].prefix == 'sky':
                # Remember the key of the target list so that it will appear in the
                # final assigned targets's list
                pd_append_column(target_list.data, '__key', k, pd.StringDtype())

                # Make room for target_idx in the target list for joining in the magnitudes later
                pd_append_column(target_list.data, '__target_idx', pd.NA, pd.Int64Dtype())

                # No cross-matching necessary here
                idx = netflow.append_sky_targets(target_list)
                target_list.data['__target_idx'] = idx

    def __cross_match_source_target_list(self, netflow, target_list):
        # Cross-match the target list with the target lists already added to netflow
        # to remove duplicates.

        if netflow.targets is None:
            return None, None, None, None
        else:
            logger.info(f'Cross-matching catalog {target_list.name} with {netflow.targets.shape[0]} already processed targets.')

            # Coordinates of the target list
            ra, dec = target_list.get_coords()
            target_coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

            # Coordinates of the targets already added to netflow
            # TODO: this is an ever increasing list so maybe cache the coordinates or even the kd-tree
            #       but the complexity of the update is probably not worth it
            ra = np.array(netflow.targets['RA'])
            dec = np.array(netflow.targets['Dec'])
            netflow_coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

            # Find the closest neighbor of each new target to the targets already added
            # idx: target list index into the netflow targets
            # sep2d: angular separation in degrees to the closest neighbor
            idx, sep2d, dist3d = match_coordinates_sky(target_coords, netflow_coords, nthneighbor=1)

            # When there are multiple matches, only keep the closest neighbors. This is
            # achieved the easiest using pandas.duplicated
            match = pd.DataFrame({'idx': idx, 'sep2d': sep2d.arcsec})

            # Filter out targets that are too far away from the closest neighbor
            mask = (match['sep2d'] <= self.__xmatch_rad * u.arcsecond)

            logger.info(f'Cross-matching resulted in {mask.sum()} matches with separation less than {self.__xmatch_rad} arc sec.')

            # Sort the good matches by idx and then separation
            match = match[mask].sort_values(by=['idx', 'sep2d'])

            # Remove the duplicated that aren't the closest match
            duplicated = match.duplicated(subset='idx', keep='first')
            match = match[~duplicated]

            logger.info(f'Removed {duplicated.sum()} duplicated matches, only keeping closest neighbors.')

            # Sort the matches by the index column
            match.sort_index(inplace=True)

            idx1 = match.index.values
            idx2 = match['idx'].values
            sep = match['sep2d'].values
            
            # Generate a mask that selects the targets that have matches
            mask = np.full(target_list.data.shape[0], False)
            mask[idx1] = True

            logger.info(f'Cross-matching resulted in {mask.sum()} unique matches with separation less than {self.__xmatch_rad} arc sec.')

            return idx1, idx2, sep, mask

if __name__ == '__main__':
    script = ImportScript()
    script.execute()