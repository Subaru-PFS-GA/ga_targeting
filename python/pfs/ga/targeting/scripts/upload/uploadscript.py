import os
import shutil
from glob import glob
import commentjson as json
import numpy as np
import pytz
from astropy.table import Table
from astropy.time import Time, TimeDelta

import ics.cobraOps
import pfs.datamodel
import pfs.utils
import pfs.instdata
import ics.cobraCharmer
from gurobipy import gurobi
from pfs.datamodel import TargetType

from ...config.netflow import NetflowConfig
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_FIELDS
from ..targetingscript import TargetingScript

from ...setup_logger import logger

class UploadScript(TargetingScript):
    """
    Command-line script to convert the netflow output into a format that can be uploaded
    to the spt_ssp_observation repository.
    """

    def __init__(self):
        super().__init__()

        self.__indir = None
        self.__outdir = None
        self.__nan_values = True
        self.__obs_run = 'march2025'
        self.__obs_wg = 'GA'
        
        self.__input_args = None            # arguments loaded from the input directory

    def _add_args(self):
        super()._add_args()

        self.add_arg('--in', type=str, required=True, help='Input directory containing the netflow output.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')

        self.add_arg('--nan-values', action='store_true', dest='nan_values', help='Use NaN values in the output files.')
        self.add_arg('--no-nan-values', action='store_false', dest='nan_values', help='Do not use NaN values in the output files.')

        self.add_arg('--obs-run', type=str, default='march2025', help='Observation run name.')
        self.add_arg('--obs-wg', type=str, default='GA', help='SSP working group abbreviation.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self.__indir = self.get_arg('in', args, self.__indir)
        self.__outdir = self.get_arg('out', args, self.__outdir)
        
    def _init_from_args(self, args):
        super()._init_from_args(args)

        self.__nan_values = self.get_arg('nan_values', args, self.__nan_values)
        self.__obs_run = self.get_arg('obs_run', args, self.__obs_run)
        self.__obs_wg = self.get_arg('obs_wg', args, self.__obs_wg)

        # Load the arguments file from the input directory and parse a few arguments
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
            self._field = M31_FIELDS[self.get_arg('m31', self.__input_args)]

        self._nvisits = self.get_arg('nvisits', self.__input_args, self._nvisits)

        if self._nvisits is not None:
            self._config.field.nvisits = self._nvisits

    def __find_last_dumpfile(self, pattern):
        """
        Find the last dump file in the input directory that matches the given pattern.
        """
        files = glob(os.path.join(self.__indir, pattern))
        files.sort()
        return files[-1] if len(files) > 0 else None

    def prepare(self):
        """
        Prepare the script for execution.
        """
        super().prepare()

        # Create the output directory
        if os.path.isdir(self.__outdir):
            raise Exception(f'Output directory already exists: {self.__outdir}')
        else:
            logger.debug(f'Creating output directory `{self.__outdir}`')
            os.makedirs(self.__outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self.__outdir, os.path.basename(self.log_file))

    def run(self):
        # Load instrument calibration data
        instrument = self._create_instrument()

        # Generate the list of pointings
        pointings = self._generate_pointings()

        # Load the final, merged assignments list
        assignments_all = self._load_assignments_all()

        # Save the list of pointings and netflow options
        self.__write_ppcList(pointings)

        # Save the list of targets
        self.__write_assignments(pointings, assignments_all)

    def _get_assignments_all_path(self):
        return os.path.join(self.__indir, f'{self._config.field.key}_assignments_all.feather')

    def __get_field_code(self, pidx, vidx):
        return f'SSP_{self.__obs_wg}_{self._config.field.key}_P{pidx:02d}_V{vidx:02d}'
    
    def __write_ppcList(self, pointings):
        """
        Write the ppcList to a file in the output directory.
        """
        
        # Give a name to each visit
        # TODO: bring out pattern as a parameter
        code = [ self.__get_field_code(pidx, vidx) for pidx in range(len(pointings)) for vidx in range(pointings[pidx].nvisits) ]

        # Pointing centers for each visit
        ra = np.array([ p.ra for p in pointings for v in range(p.nvisits) ], dtype=np.float64)
        dec = np.array([ p.dec for p in pointings for v in range(p.nvisits) ], dtype=np.float64)
        pa = np.array([ p.posang for p in pointings for v in range(p.nvisits) ], dtype=np.float32)
        resolution = np.array([ self._config.field.resolution.upper() for p in pointings for v in range(p.nvisits) ])     # TODO: bring out as parameter
        priority = np.array([ 1 for p in pointings for v in range(p.nvisits) ], dtype=np.int32)    # TODO: bring out as parameter
        
        # TODO: time zone should come from the config
        # TODO: obs_time could be different for each visit of each pointing
        #       to approximate final focal plane positions better
        hawaii_tz = pytz.timezone('US/Hawaii')
        obstime = np.array([ Time(p.obs_time.to_datetime().astimezone(hawaii_tz)).iso for p in pointings for v in range(p.nvisits) ])
        
        exptime = np.array([ int(p.exp_time.value) for p in pointings for v in range(p.nvisits) ], dtype=int)
        nframes = np.array([ 2 for p in pointings for v in range(p.nvisits) ], dtype=int)    # TODO: bring out as parameter
        
        table = Table()

        # Add package versions to the meta data

        _, _, tag = self.get_last_git_commit(pfs.datamodel)
        table.meta['datamodel'] = tag

        _, _, tag = self.get_last_git_commit(pfs.utils)
        table.meta['pfs_utils'] = tag

        _, _, tag = self.get_last_git_commit(ics.cobraCharmer)
        table.meta['ics_cobraCharmer'] = tag

        # TODO: cobraOps has no tags
        # _, _, tag = self.get_last_git_commit(ics.cobraOps)
        # table.meta['ics_cobraOps'] = tag

        # _, _, tag = self.get_last_git_commit(ets.fiberalloc)
        table.meta['ets_fiberalloc'] = None

        table.meta['gurobi'] = '{}.{}.{}'.format(*gurobi.version())

        # _, _, tag = self.get_last_git_commit(pfs.utils)
        # table.meta['Bench'] = tag

        _, _, tag = self.get_last_git_commit(pfs.instdata)
        table.meta['pfs_instdata'] = tag

        # Add runtime options to the meta data

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
            # PreSOS2Encoding, PreSOS1Encoding,
            0, 0,
            self._config.gurobi_options.threads,
        ]

        # Add columns to the table

        table['ppc_code'] = code
        table['ppc_ra'] = ra
        table['ppc_dec'] = dec
        table['ppc_pa'] = pa
        table['ppc_resolution'] = resolution
        table['ppc_priority'] = priority
        table['ppc_obstime'] = obstime
        table['ppc_exptime'] = exptime
        table['ppc_nframes'] = nframes

        fn = os.path.join(self.__outdir, f'{self.__obs_run}/targets/{self.__obs_wg}', 'ppcList.ecsv')
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        table.write(fn)

    def __get_flux_by_band(self, assignments_all):
        filter_map = {
            'g_ps1': 'g',
            'r_ps1': 'r',
            'i_ps1': 'i',
            'z_ps1': 'z',
            'y_ps1': 'y',
            'g_gaia': 'g',
            'bp_gaia': 'b',
            'rp_gaia': 'r',
        }

        bands = 'bgrizy'

        filter = { b: len(assignments_all) * [None,] for b in bands }
        psf_flux = { b: len(assignments_all) * [None,] for b in bands }
        psf_flux_err = { b: len(assignments_all) * [None,] for b in bands }

        for i, (_, ff, flux, flux_err) in enumerate(assignments_all[['filter', 'psf_flux', 'psf_flux_err']].to_records()):
            for j, f in enumerate(ff):
                b = filter_map[f]
                if filter[b][i] is None:
                    filter[b][i] = f
                    psf_flux[b][i] = flux[j]
                    psf_flux_err[b][i] = flux_err[j]

        return filter, psf_flux, psf_flux_err

    def __write_assignments(self, pointings, assignments_all):

        for target_type in [ TargetType.SCIENCE, TargetType.FLUXSTD, TargetType.SKY ]:
            for pidx in range(len(pointings)):
                for vidx in range(pointings[pidx].nvisits):
                    field_code = self.__get_field_code(pidx, vidx)

                    # Get the list of targets for this pointing and visit
                    mask = ((assignments_all['target_type'] == TargetType.SCIENCE) &
                            (assignments_all['pointing_idx'] == pidx) &
                            (assignments_all['visit_idx'] == vidx))

                    table = Table()

                    # Target list meta data

                    table.meta['ppc_code'] = field_code
                    table.meta['ppc_ra'] = pointings[pidx].ra
                    table.meta['ppc_dec'] = pointings[pidx].dec
                    table.meta['ppc_pa'] = pointings[pidx].posang

                    # Target list columns

                    obj_id = np.array(assignments_all['targetid'][mask], dtype=np.int64)
                    ra = np.array(assignments_all['RA'][mask], dtype=np.float64)
                    dec = np.array(assignments_all['Dec'][mask], dtype=np.float64)
                    pmra = np.array(assignments_all['pmra'][mask], dtype=np.float64)
                    pmdec = np.array(assignments_all['pmdec'][mask], dtype=np.float64)
                    parallax = np.array(assignments_all['parallax'][mask], dtype=np.float64)
                    epoch = np.array([f'J{e:0.2f}' for e in assignments_all['epoch'][mask] ], dtype=str)
                    tract = np.array(assignments_all['tract'][mask], dtype=np.int64)
                    patch = np.array(assignments_all['patch'][mask], dtype=str)
                    # catalog_id = 
                    target_type_id = np.array(assignments_all['target_type'][mask], dtype=np.int32)
                    input_catalog_id = np.array(assignments_all['targetid'][mask], dtype=np.int64)
                    ob_code = np.array(assignments_all['obcode'][mask], dtype=str)
                    proposal_id = np.array(assignments_all['proposalid'][mask], dtype=str)
                    priority = np.array(assignments_all['priority'][mask], dtype=np.int32)
                    effective_exptime = np.array(assignments_all['exp_time'][mask], dtype=np.float32)
                    fiberid = np.array(assignments_all['cobraid'][mask], dtype=np.int32)
                    pfi_x = np.array(assignments_all['fp_x'][mask], dtype=np.float64)
                    pfi_y = np.array(assignments_all['fp_y'][mask], dtype=np.float64)

                    # Group available fluxes into bands with filter names for each target
                    filter, psf_flux, psf_flux_err = self.__get_flux_by_band(assignments_all[mask])                    

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
                    table['effective_exptime'] = effective_exptime

                    for b in filter:
                        if np.any([ f is not None for f in filter[b] ]):
                            table[f'filter_{b}'] = filter[b]
                            table[f'psf_flux_{b}'] = psf_flux[b]
                            table[f'psf_flux_err_{b}'] = psf_flux_err[b]

                    table['fiberid'] = fiberid                      # This is in fact the cobra ID
                    table['pfi_X'] = pfi_x
                    table['pfi_Y'] = pfi_y

                    fn = os.path.join(
                        self.__outdir,
                        f'{self.__obs_run}/targets/{self.__obs_wg}/{str(target_type).lower()}',
                        f'{field_code}.ecsv')
                    os.makedirs(os.path.dirname(fn), exist_ok=True)
                    table.write(fn)

        def __copy_pfsDesign(self):
            dir = os.path.join(self.__outdir, f'{self.__obs_run}/design/{self.__obs_wg}')
            os.makedirs(dir, exist_ok=True)
            files = glob(os.path.join(self.__indir, 'pfsDesign*.fits'))
            for fn in files:
                logger.info(f'Copying {fn} to {dir}')
                # Copy the file to the output directory
                shutil.copy(fn, dir)

if __name__ == '__main__':
    script = UploadScript()
    script.execute()