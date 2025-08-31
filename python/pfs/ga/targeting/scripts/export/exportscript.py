import os
import shutil
from glob import glob
import commentjson as json
import numpy as np
import pandas as pd
import pytz
from astropy.table import Table
from astropy.time import Time, TimeDelta
import astropy.units as u

import ics.cobraOps
import pfs.datamodel
import pfs.utils
import pfs.instdata
import ics.cobraCharmer
import pfs.ga.targeting
from gurobipy import gurobi
from pfs.datamodel import TargetType, PfsDesign

from ...config.netflow import NetflowConfig
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_SECTORS
from ..targetingscript import TargetingScript

from ...setup_logger import logger

class ExportScript(TargetingScript):
    """
    Command-line script to convert the netflow output into a format that can be uploaded
    to the spt_ssp_observation repository.
    """

    def __init__(self):
        super().__init__()

        self.__nan_values = True
        self.__obs_run = '2025-03'
        self.__obs_wg = 'GA'
        self.__nframes = 4                  # Number of frame repeats (for CR removal)
        self.__input_catalog_id = None      # Overrides whatever catId is used during fiber assignment
        self.__proposal_id = None           # Overrides whatever proposalId is used during fiber assignment
        
        self.__input_args = None            # arguments loaded from the input directory

    def _add_args(self):
        super()._add_args()

        self.add_arg('--in', type=str, nargs='+', required=True, help='Input directories containing the netflow output.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')

        self.add_arg('--nan-values', action='store_true', dest='nan_values', help='Use NaN values in the output files.')
        self.add_arg('--no-nan-values', action='store_false', dest='nan_values', help='Do not use NaN values in the output files.')

        self.add_arg('--obs-run', type=str, help='Observation run name.')
        self.add_arg('--obs-wg', type=str, help='SSP working group abbreviation.')
        self.add_arg('--nframes', type=int, help='Number of frame repeats for CR removal.')
        self.add_arg('--input-catalog-id', type=int, help='Overrides whatever catId is used during fiber assignment.')
        self.add_arg('--proposal-id', type=str, help='Overrides whatever proposalId is used during fiber assignment.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')    

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self._indir = sorted(self.get_arg('in', args, self._indir))
        self._outdir = self.get_arg('out', args, self._outdir)
        
    def _init_from_args(self, args):
        super()._init_from_args(args)

        self.__nan_values = self.get_arg('nan_values', args, self.__nan_values)
        self.__obs_run = self.get_arg('obs_run', args, self.__obs_run)
        self.__obs_wg = self.get_arg('obs_wg', args, self.__obs_wg)
        self.__nframes = self.get_arg('nframes', args, self.__nframes)
        self.__input_catalog_id = self.get_arg('input_catalog_id', args, self.__input_catalog_id)
        self.__proposal_id = self.get_arg('proposal_id', args, self.__proposal_id)

        # Load the arguments file from the first input directory and parse a few arguments
        fn = self.__find_last_dumpfile('*.args')
        logger.info(f'Using input args file `{fn}`')
        with open(fn) as f:
            self.__input_args = json.loads(f.read())

        fn = self.__find_last_dumpfile('*.config')
        logger.info(f'Using input config file `{fn}`')
        self._config = NetflowConfig.from_file(fn, ignore_collisions=True, format='.json')

        if self.is_arg('dsph', self.__input_args):
            self._field = DSPH_FIELDS[self.get_arg('dsph', self.__input_args)]

        if self.is_arg('m31', self.__input_args):
            self._field = M31_SECTORS[self.get_arg('m31', self.__input_args)]

    def __find_last_dumpfile(self, pattern):
        """
        Find the last dump file in the input directory that matches the given pattern.
        """
        files = glob(os.path.join(self._indir[-1], pattern))
        files.sort()
        return files[-1] if len(files) > 0 else None

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
        # Load instrument calibration data
        # instrument = self._create_instrument()

        # Load the list of designs
        designs = self.__load_design_list()

        # Load the final, merged assignments list
        assignments = self._load_fiber_assignments_all_list()

        # Save the list of pointings and netflow options
        self.__write_ppcList(designs)

        # Save the list of targets
        self.__write_assignments(designs, assignments)

        # Copy the pfsDesign files to the output directory
        self.__copy_pfsDesigns(designs)

    def _load_fiber_assignments_all_list(self):
        # Load the list of fiber assigments from each input directory
        fiber_assignemnts_all_list = None
        for dir in self._indir:
            fn = self._get_fiber_assignments_all_path(dir)
            df = pd.read_feather(fn)

            if fiber_assignemnts_all_list is None:
                fiber_assignemnts_all_list = df
            else:
                fiber_assignemnts_all_list = pd.concat([fiber_assignemnts_all_list, df])

        return fiber_assignemnts_all_list

    def __load_design_list(self):
        # Load the list of design files from each input directory
        design_list = None
        for dir in self._indir:
            fn = self._get_design_list_path(dir)
            df = pd.read_feather(fn)

            if design_list is None:
                design_list = df
            else:
                design_list = pd.concat([design_list, df])

        return design_list

    def _get_design_list_path(self, dir):
        return os.path.join(dir, f'{self._config.field.key}_designs.feather')
    
    def __get_field_code(self, stage, pidx, vidx):
        if stage is None:
            return f'SSP_{self.__obs_wg}_{self._config.field.key}_P{pidx:02d}V{vidx:02d}'
        else:
            return f'SSP_{self.__obs_wg}_{self._config.field.key}_S{stage:01d}P{pidx:02d}V{vidx:02d}'

    def __get_ppcList_path(self):
        return os.path.join(self._outdir, f'runs/{self.__obs_run}/targets/{self.__obs_wg}', 'ppcList.ecsv')
    
    def __write_ppcList(self, designs):
        """
        Write the ppcList to a file in the output directory.
        """

        table = Table()
        
        # Give a name to each visit
        code = designs[['stage', 'pointing_idx', 'visit_idx']].apply(lambda x: self.__get_field_code(x['stage'], x['pointing_idx'], x['visit_idx']), axis=1)
        resolution = len(designs) * [ self._config.field.resolution.upper() ]

        # Convert to HST
        hawaii_tz = TimeDelta(-10 * u.hr)
        obstime = designs['obs_time'].apply(lambda t: (Time(t) + hawaii_tz).iso.replace(' ', 'T'))

        nframes = code.size * [ self.__nframes ]

        table['ppc_code'] = np.array(code, dtype=str)
        table['ppc_ra'] = np.array(designs['ra'], dtype=np.float64)
        table['ppc_dec'] = np.array(designs['dec'], dtype=np.float64)
        table['ppc_pa'] = np.array(designs['posang'], dtype=np.float32)
        table['ppc_resolution'] = np.array(resolution, dtype=str)
        table['ppc_priority'] = np.array(designs['priority'], dtype=np.int32)
        table['ppc_obstime'] = np.array(obstime, dtype=str)
        table['ppc_exptime'] = designs['exp_time']
        table['ppc_nframes'] = np.array(nframes, dtype=np.int32)
        table['ppc_pfsDesignId'] = np.array(designs['pfsDesignId'].apply(lambda id: f'0x{id:016x}'), dtype=str)

        # Add package versions to the meta data

        _, _, tag = self.get_last_git_commit(pfs.datamodel)
        table.meta['datamodel'] = tag

        _, _, tag = self.get_last_git_commit(pfs.utils)
        table.meta['pfs_utils'] = tag

        _, _, tag = self.get_last_git_commit(ics.cobraCharmer)
        table.meta['ics_cobraCharmer'] = tag

        _, _, tag = self.get_last_git_commit(ics.cobraOps)
        table.meta['ics_cobraOps'] = tag

        # _, _, tag = self.get_last_git_commit(ets.fiberalloc)
        table.meta['ets_fiberalloc'] = None

        table.meta['gurobi'] = '{}.{}.{}'.format(*gurobi.version())

        # _, _, tag = self.get_last_git_commit(pfs.utils)
        # table.meta['Bench'] = tag

        _, _, tag = self.get_last_git_commit(pfs.instdata)
        table.meta['pfs_instdata'] = tag

        hash, _, _ = self.get_last_git_commit(pfs.ga.targeting)
        table.meta['ga_targeting'] = hash

        # Add runtime options to the meta data

        table.meta['cobra safety margin'] = self._config.netflow_options.cobra_safety_margin
        table.meta['cobra maximum distance'] = self._config.netflow_options.cobra_maximum_distance
        table.meta['dot margin'] = self._config.instrument_options.black_dot_radius_margin
        table.meta['dot penalty'] = self._config.netflow_options._NetflowOptionsConfig__black_dot_penalty_str
        table.meta['gurobi_parameters'] = [
            self._config.gurobi_options.seed,
            self._config.gurobi_options.presolve,
            self._config.gurobi_options.method,
            self._config.gurobi_options.degenmoves,
            self._config.gurobi_options.heuristics,
            self._config.gurobi_options.mipfocus,
            self._config.gurobi_options.mipgap,
            0, 0,   # PreSOS2Encoding, PreSOS1Encoding,
            self._config.gurobi_options.threads,
        ]

        fn = self.__get_ppcList_path()
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        table.write(fn)

    def __get_flux_by_band(self, assignments_all):

        # TODO: move these dictionaries to some kind of settings file or config section and make
        #       them the default
        filter_bands = {
            'g_ps1': 'g',
            'r_ps1': 'r',
            'i_ps1': 'i',
            'z_ps1': 'z',
            'y_ps1': 'y',
            'g_gaia': 'g',
            'bp_gaia': 'b',
            'rp_gaia': 'r',
            'g_hsc': 'g',
            'r_hsc': 'r',
            'i_hsc': 'i',
        }

        # TODO: is this map valid for all GA fields?
        if self._field is not None:
            filter_map = self._field.get_filter_map()
        else:
            filter_map = {}

        bands = 'bgrizy'

        filter = { b: len(assignments_all) * [None,] for b in bands }
        psf_flux = { b: len(assignments_all) * [None,] for b in bands }
        psf_flux_err = { b: len(assignments_all) * [None,] for b in bands }

        for i, (_, ff, flux, flux_err) in enumerate(assignments_all[['filter', 'psf_flux', 'psf_flux_err']].to_records()):
            for j, f in enumerate(ff):
                if f in filter_bands:
                    b = filter_bands[f]
                    if b is not None and filter[b][i] is None:
                        filter[b][i] = filter_map[f] if f in filter_map else f        ####
                        psf_flux[b][i] = flux[j]
                        psf_flux_err[b][i] = flux_err[j]

        return filter, psf_flux, psf_flux_err

    def __get_targets_path(self, target_type, field_code):
        return os.path.join(
            self._outdir,
            f'runs/{self.__obs_run}/targets/{self.__obs_wg}/{str(target_type).lower()}',
            f'{field_code}.ecsv')

    def __write_assignments(self, designs, assignments):

        for target_type in [ TargetType.SCIENCE, TargetType.FLUXSTD, TargetType.SKY ]:
            for i, (stage, pidx, vidx) in designs[['stage', 'pointing_idx', 'visit_idx']].iterrows():
                field_code = self.__get_field_code(stage, pidx, vidx)

                table = Table()
                
                # Target list meta data

                mask = ((designs['pointing_idx'] == pidx) &
                        (designs['visit_idx'] == vidx))
                
                if stage is not None:
                    mask &= (designs['stage'] == stage)

                table.meta['ppc_code'] = field_code
                table.meta['ppc_ra'] = designs[mask]['ra'].item()
                table.meta['ppc_dec'] = designs[mask]['dec'].item()
                table.meta['ppc_pa'] = designs[mask]['posang'].item()

                # Get the list of targets for this pointing and visit
                mask = ((assignments['target_type'] == target_type) &
                        (assignments['pointing_idx'] == pidx) &
                        (assignments['visit_idx'] == vidx))

                if stage is not None:
                    mask &= (assignments['stage'] == stage)

                assert mask.sum() > 0

                # Target list columns
                if target_type == TargetType.SCIENCE:
                    id_prefix = self._get_id_prefix()
                    obj_id = (id_prefix | np.array(assignments['__target_idx'][mask], dtype=np.int64))
                else:
                    obj_id = np.array(assignments['targetid'][mask], dtype=np.int64)

                ra = np.array(assignments['RA'][mask], dtype=np.float64)
                dec = np.array(assignments['Dec'][mask], dtype=np.float64)
                pmra = np.array(assignments['pmra'][mask].fillna(0.0), dtype=np.float64)
                pmdec = np.array(assignments['pmdec'][mask].fillna(0.0), dtype=np.float64)
                parallax = np.array(assignments['parallax'][mask].fillna(1.0e-5), dtype=np.float64)
                epoch = np.array([f'J{e:0.2f}' for e in assignments['epoch'][mask] ], dtype=str)
                tract = np.array(assignments['tract'][mask].fillna(-1), dtype=np.int64)
                patch = np.array(assignments['patch'][mask].fillna('-1,-1'), dtype=str)
                # catalog_id = 
                target_type_id = np.array(assignments['target_type'][mask], dtype=np.int32)

                if target_type == TargetType.SCIENCE and self.__input_catalog_id is not None:
                    input_catalog_id = np.full(mask.sum(), self.__input_catalog_id, dtype=np.int64)
                else:
                    input_catalog_id = np.array(assignments['catid'][mask], dtype=np.int64)

                ob_code = np.array(assignments['obcode'][mask], dtype=str)
                
                if self.__proposal_id is not None:
                    proposal_id = np.array(mask.sum() * [self.__proposal_id], dtype=str)
                else:
                    proposal_id = np.array(assignments['proposalid'][mask], dtype=str)
                
                if target_type == TargetType.SCIENCE:
                    priority = np.array(assignments['priority'][mask], dtype=np.int32)
                else:
                    priority = np.full(mask.sum(), -1, dtype=np.int32)

                if target_type == TargetType.FLUXSTD:
                    if 'prob_f_star' in assignments:
                        prob_f_star = np.array(assignments['prob_f_star'][mask], dtype=float)
                    else:
                        prob_f_star = np.full(mask.sum(), 1, dtype=float)
                else:
                    prob_f_star = np.full(mask.sum(), -1, dtype=float)

                effective_exptime = np.array(assignments['exp_time'][mask].fillna(0.0), dtype=np.float32)
                cobraid = np.array(assignments['cobraid'][mask], dtype=np.int32)
                pfi_x = np.array(assignments['fp_x'][mask], dtype=np.float64)
                pfi_y = np.array(assignments['fp_y'][mask], dtype=np.float64)

                # Group available fluxes into bands with filter names for each target
                filter, psf_flux, psf_flux_err = self.__get_flux_by_band(assignments[mask])                    

                # Put together the table

                table['obj_id'] = obj_id
                table['ra'] = ra
                table['dec'] = dec
                table['pmra'] = pmra
                table['pmdec'] = pmdec
                table['parallax'] = parallax
                table['epoch'] = epoch
                table['tract'] = tract
                table['patch'] = patch
                # table['catalog_id'] = catalog_id
                table['target_type_id'] = target_type_id
                table['input_catalog_id'] = input_catalog_id
                table['ob_code'] = ob_code
                table['proposal_id'] = proposal_id
                table['priority'] = priority
                table['prob_f_star'] = prob_f_star
                table['effective_exptime'] = effective_exptime

                for b in filter:
                    if np.any([ f is not None for f in filter[b] ]):
                        table[f'filter_{b}'] = np.array([f if f is not None else '' for f in filter[b]], dtype=str)
                        table[f'psf_flux_{b}'] = np.array([f if f is not None else np.nan for f in psf_flux[b]], dtype=float)
                        table[f'psf_flux_error_{b}'] = np.array([f if f is not None else np.nan for f in psf_flux_err[b]], dtype=float)
                    else:
                        table[f'filter_{b}'] = np.array(mask.sum() * [''], dtype=str)
                        table[f'psf_flux_{b}'] = np.full(mask.sum(), np.nan, dtype=float)
                        table[f'psf_flux_error_{b}'] = np.full(mask.sum(), np.nan, dtype=float)

                table['cobraId'] = cobraid
                table['pfi_X'] = pfi_x
                table['pfi_Y'] = pfi_y

                fn = self.__get_targets_path(target_type, field_code)
                os.makedirs(os.path.dirname(fn), exist_ok=True)
                table.write(fn)

    def __get_pfsDesign_path(self):
        return os.path.join(self._outdir, f'runs/{self.__obs_run}/pfs_designs/{self.__obs_wg}')

    def __copy_pfsDesigns(self, designs):

        # Create directory for output files
        outdir = self.__get_pfsDesign_path()
        os.makedirs(outdir, exist_ok=True)

        for i, (stage, pidx, vidx, pfsDesignId) in designs[['stage', 'pointing_idx', 'visit_idx', 'pfsDesignId']].iterrows():
            # Find the design in one of the input directories
            fn = PfsDesign.fileNameFormat % (pfsDesignId)
            dir = None
            for indir in self._indir:
                if os.path.isfile(os.path.join(indir, fn)):
                    dir = indir
                    break

            if dir is None:
                raise FileNotFoundError(f'Design file {fn} not found.')

            # Load the design file and override name
            d = PfsDesign.read(pfsDesignId, dirName=dir)
            d.designName = self.__get_field_code(stage, pidx, vidx)

            # Save the design to the output directory
            logger.info(f'Copying {fn} to {outdir}')
            d.write(dirName=outdir)

if __name__ == '__main__':
    script = ExportScript()
    script.execute()