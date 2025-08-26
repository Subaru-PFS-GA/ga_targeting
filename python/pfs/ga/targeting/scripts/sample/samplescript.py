import os
import numpy as np
from scipy.special import logsumexp

from pfs.ga.isochrones.isogrid import IsoGrid
import pfs.ga.targeting
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_SECTORS
from ...config.sample import SampleConfig
from ..script import Script
from ...io import Hdf5SimulationReader
from ...instrument import SubaruHSC
from ... import ProbabilityMap
from ...selection import ProbabilityCut, ProbabilitySampling, MagnitudeSelection, ColorSelection, LinearSelection
from ...io import GaiaReader, ObservationSerializer

from ...setup_logger import logger

class SampleScript(Script):
    """
    Script to apply selections and assign priorities to the targets based on HSC data.
    """

    # TODO: merge repeated functions with PMapScript

    def __init__(self):
        super().__init__()

        self._field = None
        self._config = None

        self.__outdir = None
        self.__gaiadir = None
        self.__skip_notebooks = False

    def __get_config(self):
        return self._config
    
    config = property(__get_config)

    def _add_args(self):
        super()._add_args()

        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')
        self.add_arg('--gaiadir', type=str, required=False, help='Path to the directory containing the Gaia HDF5 file.')

        self.add_arg('--dsph', type=str, choices=DSPH_FIELDS.keys(), help='Name of a predefined dSph target.')
        self.add_arg('--m31', type=str, choices=M31_SECTORS, help='Name of a predefined M31 field.')

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self.__outdir = self.get_arg('out', args, self.__outdir)
        if self.is_arg('gaiadir'):
            self.__gaiadir = self.get_arg('gaiadir', args, self.__gaiadir)
        else:
            self.__gaiadir = os.path.join(self.__outdir, '../../gaia')

    def _init_from_args(self, args):
        super()._init_from_args(args)

        if self.is_arg('dsph', args):
            self._field = DSPH_FIELDS[self.get_arg('dsph', args)]

        if self.is_arg('m31', args):
            self._field = M31_SECTORS[self.get_arg('m31', args)]

        # If a field is specified, load its default configuration      
        if self._field is not None:
            self._config = self._field.get_sample_config()
        else:
            self._config = SampleConfig.default()

        # Load the configuration template files and merge with the default config
        config_files = self.get_arg('config', args)
        self._config.load(config_files, ignore_collisions=True)
        logger.info(f'Loaded {len(config_files)} config files from {config_files}.')

        self.__skip_notebooks = self.get_arg('skip_notebooks', args, self.__skip_notebooks)

    def prepare(self):
        super().prepare()

        # Create the output directory
        if os.path.isdir(self.__outdir):
            raise Exception(f'Output directory already exists: {self.__outdir}')
        else:
            logger.debug(f'Creating output directory `{self.__outdir}`')
            os.makedirs(self.__outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self.__outdir, os.path.basename(self.log_file))

    def __get_output_config_path(self):
        command = self.get_command_name()
        path = os.path.join(self.__outdir, f'{command}_{self.timestamp}.config')
        return path

    def _dump_settings(self):
        super()._dump_settings()

        # Save the active configuration to the output directory
        path = self.__get_output_config_path()
        self._config.save(path)

        logger.debug(f'Configuration saved to `{os.path.abspath(path)}`.')

    def run(self):
        hsc = self.__load_observations()
        pmap = self.__load_pmap()
        iso = self.__load_isochrones()

        if self._config.gaia_crossmatch:
            if self.__gaia_exists():
                gaia = self.__load_gaia()
            else:
                gaia = self.__query_gaia()
                self.__save_gaia(gaia)

            hsc_gaia_idx, hsc_gaia_mask, hsc_gaia_sep = self.__crossmatch_hsc_gaia(hsc, gaia)
            gaia_hsc_idx, gaia_hsc_mask, gaia_hsc_sep = self.__crossmatch_gaia_hsc(hsc, gaia)
            self.__merge_proper_motions(hsc, gaia, hsc_gaia_idx, hsc_gaia_mask)
        else:
            gaia = None

        # Define the cut based on membership probability
        probcut = ProbabilityCut(pmap, 1, self._config.lp_member_limit)

        # Apply the selection and assing probabilites
        selection, mask = self.__assign_probabilities(hsc, pmap, probcut)

        self.__assign_priorities(hsc, isogrid=iso)
        self.__save_target_list(hsc)

        # Execute the evaluation notebooks
        if not self.__skip_notebooks:
            for notebook in ['sample']:
                logger.info(f'Executing evaluation notebook `{notebook}`...')
                notebook_path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), f'scripts/sample/notebooks/{notebook}.ipynb')
                parameters = {
                    'DEBUG': False,
                    'CONFIG_FILE': self.__get_output_config_path(),
                    'OUTPUT_PATH': self.__outdir,
                }
                self._execute_notebook(notebook_path, parameters, self.__outdir)

    def __load_observations(self):
        reader = self._field.get_text_observation_reader()
        obs = reader.read(os.path.expandvars(self._config.obs_path))
        return obs

    def __load_pmap(self):
        cmd = self._field.get_cmd()
        pmap = ProbabilityMap(cmd.axes)
        pmap.load(os.path.expandvars(os.path.join(self._config.pmap_path, 'pmap.h5')))
        return pmap

    def __load_isochrones(self):
        fn = os.path.expandvars(os.path.join(self._config.isochrones_path, 'isochrones.h5'))
        iso = IsoGrid()
        iso.load(fn)

        logger.info(f'Loaded isochrones from `{fn}`.')

        return iso

    def __query_gaia(self):
        logger.info('Querying Gaia DR3...')
        reader = GaiaReader()
        gaia = reader.cone_search((self._field.pos.ra, self._field.pos.dec), self._field.rad)

        logger.info(f'Found {gaia.shape[0]} stars in the GAIA database.')

        return gaia

    def __get_gaia_file_path(self):
        return os.path.expandvars(os.path.join(self.__gaiadir, f'gaia_{self._field.ID}.feather'))

    def __gaia_exists(self):
        fn = self.__get_gaia_file_path()
        return os.path.isfile(fn)

    def __load_gaia(self):
        fn = self.__get_gaia_file_path()
        r = ObservationSerializer()
        gaia = r.read(fn)
        return gaia

    def __save_gaia(self, gaia):
        fn = self.__get_gaia_file_path()
        w = ObservationSerializer()
        w.write(gaia, fn)

    def __crossmatch_hsc_gaia(self, hsc, gaia):
        """
        Cross-match the observations with GAIA and use the proper motions from the
        closest neighbors
        """

        hsc_gaia_idx, separation = hsc.cross_match(gaia)
        logger.info(f'GAIA->HSC source median separation: {np.median(separation.arcsec)} arcsec.')
        hsc_gaia_mask = (separation.arcsec < self._config.gaia_crossmatch_radius)
        logger.info(f'HSC targets with matching GAIA targets: {hsc_gaia_mask.sum()}')

        return hsc_gaia_idx, hsc_gaia_mask, separation

    def __crossmatch_gaia_hsc(self, hsc, gaia):
        gaia_hsc_idx, separation = gaia.cross_match(hsc)
        logger.info(f'HSC->GAIA source median separation: {np.median(separation.arcsec)} arcsec.')
        gaia_hsc_mask = (separation.arcsec < self._config.gaia_crossmatch_radius)
        logger.info(f'GAIA targets with matching HSC targets: {gaia_hsc_mask.sum()}')

        return gaia_hsc_idx, gaia_hsc_mask, separation

    def __merge_proper_motions(self, hsc, gaia, hsc_gaia_idx, hsc_gaia_mask):
        """
        Merge the proper motions from GAIA and HSC
        """

        logger.info('Merging proper motions from GAIA into HSC...')
        columns = ['parallax', 'pm', 'pmdec', 'pmra', 'err_parallax', 'err_pmdec', 'err_pmra']
        hsc.merge(gaia, hsc_gaia_idx, columns=columns, mask=hsc_gaia_mask)

    def __assign_probabilities(self, obs, pmap, probcut):
        """
        Assign probabilities to the targets based on the probability map
        """

        logger.info('Applying selection...')
        selection = self._field.get_selection_mask(obs, nb=self._config.cut_nb, probcut=probcut)
        logger.info(f'HSC target count: {obs.data.shape}, selected target count: {selection.sum()}.')

        logger.info('Assigning probabilities...')
        self._field.assign_probabilities(obs, pmap, mask=selection)
        logger.info(f"HSC targets with assigned probabilities: {(~np.isnan(obs.data['p_member'])).sum()}.")

        mask = selection & ~np.isnan(obs.data['p_member'])

        return selection, mask

    def __assign_priorities(self, obs, mask=None, isogrid=None):
        """
        Assign priorities to the targets based on the probabilities
        """

        logger.info('Assigning priorities...')
        self._field.assign_priorities(obs, mask=mask, isogrid=isogrid)
        logger.info(f"HSC targets with assigned priorities: {(0 <= obs.data['priority']).sum()}.")
        logger.info(f"Unique priority classes: {np.sort(obs.data['priority'].unique())}.")
        logger.info(f"Unique exposure times: {np.sort(obs.data['exp_time'].unique())}.")

    def __get_target_list_path(self):
        return os.path.expandvars(os.path.join(self.__outdir, f'hsc_{self._field.ID}_priorities.feather'))

    def __save_target_list(self, hsc):
        """
        Save the target list to the output directory
        """

        # Only save targets with assigned priorities and exposure times
        mask = (hsc.data['priority'] >= 0)

        fn = self.__get_target_list_path()
        s = ObservationSerializer()
        s.write(hsc, fn, filter=mask)
        logger.info(f'Saved target list with priorities to {fn}.')

if __name__ == '__main__':
    script = SampleScript()
    script.execute()