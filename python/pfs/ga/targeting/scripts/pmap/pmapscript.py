import os
import numpy as np

import pfs.ga.targeting
from ...targets.dsph import GALAXIES as DSPH_FIELDS
from ...targets.m31 import M31_SECTORS
from ...config.pmap import PMapConfig
from ..script import Script
from ...io import Hdf5SimulationReader
from ...instrument import SubaruHSC
from ... import ProbabilityMap

from ...setup_logger import logger

class PMapScript(Script):
    """
    Script to generate the probability map from a simulated CMD
    """

    # TODO: merge repeated functions with SampleScript

    def __init__(self):
        super().__init__()

        self._field = None
        self._use_p_stars = False
        self._config = None

        self.__outdir = None
        self.__skip_notebooks = False

    def __get_config(self):
        return self._config
    
    config = property(__get_config)

    def _add_args(self):
        super()._add_args()

        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')

        self.add_arg('--dsph', type=str, choices=DSPH_FIELDS.keys(), help='Name of a predefined dSph target.')
        self.add_arg('--m31', type=str, choices=M31_SECTORS, help='Name of a predefined M31 field.')

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self.__outdir = self.get_arg('out', args, self.__outdir)

    def _init_from_args(self, args):
        super()._init_from_args(args)

        if self.is_arg('dsph', args):
            self._field = DSPH_FIELDS[self.get_arg('dsph', args)]

        if self.is_arg('m31', args):
            self._field = M31_SECTORS[self.get_arg('m31', args)]

        # If a field is specified, load its default configuration      
        if self._field is not None:
            self._config = self._field.get_pmap_config()
        else:
            self._config = PMapConfig.default()

        # # Load the configuration template files and merge with the default config
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
        sim = self.__load_simulation()
        pmap = self.__create_pmap(sim)
        self.__save_pmap(pmap)

        # Execute the evaluation notebooks
        if not self.__skip_notebooks:
            for notebook in ['pmap']:
                logger.info(f'Executing evaluation notebook `{notebook}`...')
                notebook_path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), f'scripts/pmap/notebooks/{notebook}.ipynb')
                parameters = {
                    'DEBUG': False,
                    'CONFIG_FILE': self.__get_output_config_path(),
                    'OUTPUT_PATH': self.__outdir,
                    'BINARIES': self.is_arg('m31'),
                }
                self._execute_notebook(notebook_path, parameters, self.__outdir)

    def __load_simulation(self):
        # Load the simulation from the input directory
        fn = os.path.join(self._config.sim_path, 'sample.h5')

        r = Hdf5SimulationReader()
        cm = {}

        # Include this for loop for simulations made with older isochrone tables that have different column names
        for prefix in ['', 'obs_', 'err_', 'flux_', 'obs_flux_', 'err_flux_', 'counts_', 'obs_counts_', 'err_counts_']:
            cm[prefix + 'hsc_g2'] = prefix + 'hsc_g'
            cm[prefix + 'hsc_r2'] = prefix + 'hsc_r'

        r.column_mapping = cm
        r.append_photometry(SubaruHSC.photometry())
        logger.info(f'Loading simulation from {fn}')
        sim = r.read(fn)

        return sim

    def __apply_cuts(self, sim):
        mask = self._field.get_selection_mask(sim, observed=True, nb=self._config.cut_nb, blue=self._config.keep_blue)
        return mask

    def __update_weights(self, sim):
        if self._config.population_weights is None:
            # Count the objects in each population instead of using the
            # input weights.
            # TODO: use the mask?
            w1 = np.bincount(sim.data['g']) / sim.data['g'].shape
        else:
            # Take the new weights from the config file
            w1 = np.array(self._config.population_weights)

        # Generate new random assignments
        g1 = np.random.choice(np.arange(w1.size, dtype=int), sim.data['g'].size, p=w1)

        # Verify the new assignments
        w2 = np.bincount(g1) / g1.shape

        return w1, g1

    def __create_pmap(self, sim):
        cmd = self._field.get_cmd()
        mask = self.__apply_cuts(sim)
        w1, g1 = self.__update_weights(sim)

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(
            sim,
            bins = self._config.bins,
            extents = self._config.extents,
            merge_list = self._config.merge_list,
            population_weights = w1,
            observed = True,
            use_p_stars = self._config.use_p_stars,
            mask = mask)
        pmap.maximum_filter()

        return pmap

    def __get_pmap_path(self):
        fn = 'pmap'

        # if self._config.cut_nb:
        #     fn += '_nb'

        # if self._config.keep_blue:
        #     fn += '_blue'

        fn += '.h5'

        return os.path.join(self.__outdir, fn)

    def __save_pmap(self, pmap):
        fn = self.__get_pmap_path()
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        pmap.save(fn)
        logger.info(f'Saved probability map to {fn}')

if __name__ == '__main__':
    script = PMapScript()
    script.execute()