import os
import commentjson as json
from datetime import datetime, timedelta, tzinfo
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

class NetflowScript(TargetingScript):
    """
    Command-line script to execute netflow optimazion and generate
    the target list and design files.

    This script loads the configuration file and reads the target lists. At least
    one science target list must be provided, in addition to a mandatory calibration
    target list and a list of sky positions.

    The target lists must already contain the necessary columns for the netflow
    optimization. The column names are case sensitive.

        * Magnitudes and/or fluxes and their errors for each filter, measured for
            the PSF, fiber aperture and/or total. The column names are analyzed in
            this particular order, if any of the columns is not found, it will be
            calculated from the available columns. For example, if the fiber flux is
            missing, but the PSF flux is available, the PSF flux will be copied.
            If the flux is not available, but the magnitude is available, the flux
            will be calculated from the magnitude. AB magnitudes and fluxes in nJy are
            assumed everywhere. The magnitude and flux names, as well as the filter
            names can be listed two different ways, see below.
        * epoch: epoch of the coordinates, default is J2016.0. It has no effect if
            the proper motion components are not provided. Note that this is not the
            equinox of the coordinate system.
        * catid: int32, catalog id, default is -1 for calibration targets and sky targets.
            It must have a unique value in the database.
        * proposalid: string, proposal id, default is empty string.
        * tract: int32, tract number, computed from the coordinates.
        * patch: string, patch number, computed from the coordinates.
        * targetid: int64, must be unique for each target in the target list and across all
            target lists. Can be -1 for sky positions.        
        * obcode: string, must be unique for each target in the target list and across all
            target lists and must be unique across all visits.
        * RA, Dec, float64, target coordinates in degrees.
        * pmra, pmdec, parallax: float64, proper motion components in mas/yr and parallax in mas.
            pmra is assumed to be multiplied by cos(Dec).
        * prefix: string, prefix for the target type, must be one of `sci`, `cal`, `sky`.
        * priority: int32, priority of the target, only for science targets


    Magnitudes/fluxes and filter names can be provided in two different ways:
        1. If the configuration entry `filters` is provided, the filter names are listed in
            the configuration file as a dictionary, where the keys are the filter names.
            The columns corresponding to the various fluxes and magnitudes must be listed
            individually in the configuration.
        2. If the configuration entry `bands` is provided, the band names are listed in
            the configuration file as a dictionary, where the keys are the band names.
            The columns corresponding to the various fluxes and magnitudes must be listed
            individually in the configuration for each band. In addition, a column must be
            defined that contains the filter name for each band. This method is typically
            used for calibration target lists where the targets might come from different
            photometric catalogs.

    The pointings are defined in the configuration file but the number of visits, the
    observation time and the exposure time can be overridden with the command-line arguments
    `--nvisits`, `--obs-time` and `--exp-time`. The default is always what is in the
    configuration file.
    """

    def __init__(self):
        super().__init__()

        self.__outdir = None
        self.__resume = False
        self.__exp_time = None
        self.__obs_time = None
        self.__time_limit = None
        self.__xmatch_rad = 2           # arcseconds
        self.__skip_notebooks = False

    def _add_args(self):
        super()._add_args()

        self.add_arg('--dsph', type=str, choices=DSPH_FIELDS.keys(), help='Name of a predefined dSph target.')
        self.add_arg('--m31', type=str, choices=M31_FIELDS, help='Name of a predefined M31 field.')

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')
        self.add_arg('--resume', action='store_true', help='Resume from a previous run.')
        self.add_arg('--nvisits', type=int, help='Number of visits for each pointing.')
        self.add_arg('--exp-time', type=float, help='Exposure time per visit, in seconds.')
        self.add_arg('--obs-time', type=str, help='Observation time in ISO format in UTC.')
        self.add_arg('--time-limit', type=int, help='Time limit for the optimization in seconds.')
        
        self.add_arg('--xmatch-rad', type=float, help='Cross-match radius in arcseconds.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self.__outdir = self.get_arg('out', args, self.__outdir)
        self.__resume = self.get_arg('resume', args, self.__resume)

    def _init_from_args(self, args):
        super()._init_from_args(args)

        if self.is_arg('dsph', args):
            self._field = DSPH_FIELDS[self.get_arg('dsph', args)]

        if self.is_arg('m31', args):
            self._field = M31_FIELDS[self.get_arg('m31', args)]

        # If a field is specified, load its default configuration      
        if self._field is not None:
            self._config = self._field.get_netflow_config()
        else:
            self._config = NetflowConfig.default()

        # Load the configuration template files and merge with the default config
        config_files = self.get_arg('config', args)
        self._config.load(config_files, ignore_collisions=True)
        logger.info(f'Loaded {len(config_files)} config files from {config_files}.')
        logger.info(f'Found {len(self._config.targets)} target lists in the configuration.')

        self._nvisits = self.get_arg('nvisits', args, self._nvisits)
        self.__exp_time = self.get_arg('exp_time', args, self.__exp_time)
        self.__obs_time = self.get_arg('obs_time', args, self.__obs_time)
        self.__time_limit = self.get_arg('time_limit', args, self.__time_limit)
        self.__xmatch_rad = self.get_arg('xmatch_rad', args, self.__xmatch_rad)
        self.__skip_notebooks = self.get_arg('skip_notebooks', args, self.__skip_notebooks)

        # Override the configuration with the command-line arguments
        if self._nvisits is not None:
            self._config.field.nvisits = self._nvisits

        if self.__exp_time is not None:
            self._config.field.exp_time = self.__exp_time

        if self.__obs_time is not None:
            self.__obs_time = datetime.fromisoformat(self.__obs_time)
            self._config.field.obs_time = self.__obs_time

        if self.__time_limit is not None:
            self._config.gurobi_options.timelimit = self.__time_limit

    def prepare(self):
        super().prepare()

        # Create the output directory
        if not self.__resume and os.path.isdir(self.__outdir):
            raise Exception(f'Output directory already exists: {self.__outdir}')
        elif self.__resume and not os.path.isdir(self.__outdir):
            raise Exception(f'Output directory does not exist: {self.__outdir}, cannot resume.')
        elif not self.__resume:
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

        # Load instrument calibration data
        instrument = self._create_instrument()

        # Generate the list of pointings
        pointings = self._generate_pointings()

        # Create the netflow object
        netflow = Netflow(
            f'{self._config.field.name}',
            instrument,
            pointings,
            workdir = self.__outdir,
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
        self.__append_source_target_lists(netflow, target_lists)
        self.__validate_targets(netflow)

        # Run the netflow optimization                
        self.__run_netflow(netflow)

        # We do not implement resume logic beyond this point because it's not worth it time-wise.

        # Verify collisions and unassign cobras with lower priority targets
        # TODO: this would check trajectory collisions but this test is not required
        #       when the cobras can move radially
        # self.__unassign_colliding_cobras(netflow)
        
        # Extract the assignments from the netflow solution
        # The column target_idx is a unique index into the netflow target cache
        assignments = self.__extract_assignments(netflow)
        self.__save_assignments(assignments)

        # Generate a summary of the target assignments
        # The index of the data frame is a unique index into the netflow target cache
        summary = netflow.get_target_assignment_summary()
        self.__save_summary(summary)

        # TODO: plot assignment statistics

        # Join the targets lists with the assignments to append the fluxes and
        # filter names as columns of lists, as  required for the design files
        assigned_target_lists = self.__append_flux_filter_lists(assignments, target_lists)

        # Concatenate the assigned target lists into a single DataFrame with only
        # the columns that are needed for the design files
        assigned_targets_all = self.__merge_assigned_target_lists(assigned_target_lists)

        # TODO: this is the point to postfix obcode with the visit id

        # Join the assignments with the merged target list to append the fluxes and
        # filter names, ob_code, etc.
        assignments_all = self.__join_assignments(assignments, assigned_targets_all)
        self._save_assignments_all(assignments_all)

        # Generate the designs and save them to the output directory
        designs = self.__create_designs(netflow, assignments_all)

        # Verify the final target lists
        self.__verify_designs(designs)

        # # Verify each design
        # # TODO: move this after any other step and throw exception on
        # #       verification failure, such az endpoint collisions or 
        # #       elbow collisions. Only report trajectory collisions.
        # logger.info('Verifying solution.') 
        # netflow.verify()

        self.__save_designs(designs)
        self._save_design_list(netflow, designs)

        # Execute the evaluation notebooks
        if not self.__skip_notebooks:
            for notebook in ['targets', 'calibration', 'assignments', 'cobra_groups', 'design']:
                logger.info(f'Executing evaluation notebook `{notebook}`...')
                notebook_path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), f'scripts/netflow/notebooks/{notebook}.ipynb')
                parameters = {
                    'DEBUG': False,
                    'CONFIG_FILE': self.__get_output_config_path(),
                    'OUTPUT_PATH': self.__outdir,
                }
                self._execute_notebook(notebook_path, parameters, self.__outdir)
    
    def load_source_target_lists(self):
        target_lists = {}
        for k, target_list_config in self._config.targets.items():
            # Check if the preprocessed target lists are available

            if self.__outdir is not None:
                path = self.__get_preprocessed_target_list_path(k)
            else:
                path = None

            loaded = False
            if self.__resume and self.__outdir is not None and os.path.isfile(path):
                target_lists[k] = self.__load_preprocessed_target_list(k, path)
                loaded = True

            if not loaded:
                target_list = self.load_target_list(k, target_list_config)
                target_list.name = k
                
                self.__calculate_target_list_flux(k, target_list)
                self.__append_target_list_extra_columns(k, target_list)
                self.__validate_target_list(k, target_list_config, target_list)
                
                target_lists[k] = target_list

                if self.__outdir is not None:
                    self.__save_preprocessed_target_list(k, target_lists[k], path)

        return target_lists

    @staticmethod
    def load_target_list(key, target_list_config, fn=None):
        if target_list_config.reader is not None:
            # TODO: implement custom reader function
            raise NotImplementedError()
        
        fn = fn if fn is not None else os.path.expandvars(target_list_config.path)
        
        logger.info(f'Reading target list from `{fn}`.')

        reader = ObservationSerializer(
            catalog_name=key,
            columns = target_list_config.columns,
            column_map = target_list_config.column_map,
            data_types = target_list_config.data_types,
            value_map = target_list_config.value_map,
            index = target_list_config.index,
            filters = target_list_config.photometry.filters if target_list_config.photometry is not None else None,
            bands = target_list_config.photometry.bands if target_list_config.photometry is not None else None,
            limits = target_list_config.photometry.limits if target_list_config.photometry is not None else None,
            mask = target_list_config.mask,
            kwargs = target_list_config.reader_args,
        )
        catalog = reader.read(fn)
        catalog.name = key
        
        logger.info(f'Read {catalog.shape[0]} targets from target list.')

        return catalog
    
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
        set_extra_column_constant('epoch', 'string', 'J2016.0', self._config.targets[key].epoch)
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

    def __validate_targets(self, netflow: Netflow):
        # TODO: Make sure cobra location/instrument group constraints are satisfiable
        
        # TODO: Move these under netflow and run on the cached target list
        
        netflow.validate_science_targets()
        netflow.validate_fluxstd_targets()
        netflow.validate_sky_targets()

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

        target_epoch = self._config.netflow_options.epoch
        target_epoch = normalize_epoch(target_epoch)
        
        logger.info(f'Converting target list `{target_list.name}` to ICRS and epoch J{target_epoch:0.2f}.')
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
    
    def __get_preprocessed_target_list_path(self, key):
        return os.path.join(self.__outdir, f'{self._config.field.key}_targets_{key}.feather')

    def __save_preprocessed_target_list(self, key, target_list, path):
        logger.info(f'Saving preprocessed target list `{key}` to `{path}`.')
        s = ObservationSerializer()
        s.write(target_list, path)

    def __load_preprocessed_target_list(self, key, path):
        logger.info(f'Loading preprocessed target list `{key}` list from `{path}`.')
        s = ObservationSerializer()
        target_list = s.read(path)
        target_list.name = key
        return target_list
    
    def __run_netflow(self, netflow):
        logger.info('Building netflow model.')
        netflow.build(resume=self.__resume, save=True)

        logger.info('Solving netflow model.')
        netflow.solve()

        # TODO: report netflow results or turn on logging

        logger.info('Saving netflow solution.')
        netflow.save_solution()

    def __unassign_colliding_cobras(self, netflow):
        logger.info('Resolving cobra collisions by unassigning low priority targets.')
        netflow.simulate_collisions()
        netflow.unassign_colliding_cobras()

        # Rerun the simulation to see if there's any remaining collisions
        netflow.simulate_collisions()
        netflow.simulate_collisions()

    def __extract_assignments(self, netflow):
        """
        Extract the assignments from the netflow problem as a DataFrame with all
        necessary information to generate the design files.
        """
        
        assignments = netflow.get_fiber_assignments(
            include_target_columns=True,
            include_unassigned_fibers=True,
            include_engineering_fibers=True)

        return assignments
    
    def __get_assignments_path(self):
        return os.path.join(self.__outdir, f'{self._config.field.key}_assignments.feather')
    
    def __save_assignments(self, assignments):
        path = self.__get_assignments_path()
        logger.info(f'Saving fiber assignments to `{path}`.')
        assignments.to_feather(path)

    def __load_assignments(self):
        path = self.__get_assignments_path()
        logger.info(f'Loading fiber assignments from `{path}`.')
        return pd.read_feather(path)
    
    def __get_summary_path(self):
        return os.path.join(self.__outdir, f'{self._config.field.key}_summary.feather')
    
    def __save_summary(self, summary):
        path = self.__get_summary_path()
        logger.info(f'Saving target assignment summary to `{path}`.')
        summary.to_feather(path)
    
    def _get_assignments_all_path(self):
        return os.path.join(self.__outdir, f'{self._config.field.key}_assignments_all.feather')

    def __append_flux_filter_lists(self, assignments, target_lists):
        """
        Join the assignments with the target lists to append the fluxes that will
        be included in the design files.
        """

        # Indexes of all assigned targets
        # Make it unique in case a target is assigned multiple times (multiple visits)
        assigned_target_idx = pd.DataFrame({'__target_idx': assignments[assignments['target_idx'] != -1]['target_idx'].unique()})

        assigned_target_lists = {}
        for k, config in self._config.targets.items():
            # Only include targets that can be assigned to fibers, ie. no guide stars
            if config.prefix in ['sci', 'cal', 'sky']:
                target_list = target_lists[k]

                # Filter down the target list to only include targets that were assigned
                assigned = target_list.data.set_index('__target_idx').join(assigned_target_idx.set_index('__target_idx'), how='inner')
                assigned.reset_index(inplace=True)

                # Filter names are either defined in the config or can be extracted from the columns
                # If the config entry is a dictionary, we assume that the keys are the filter names.
                # If it is a list, we assume that the filter names are stored as values in a column.

                photometry = self._config.targets[k].photometry

                if photometry is not None and photometry.filters is not None:
                    # Convert the flux columns into a single column of list objects
                    # Refer to canonical column names in __calculate_target_list_flux_filters

                    def flux_to_list(row, prefix, postfix):    
                        return [ row[f'{b}_{prefix}flux{postfix}'] for b in photometry.filters.keys() ]
                    
                    def filter_to_list(row):
                        return list(photometry.filters.keys())
                    
                    for prefix in [ 'psf_', 'fiber_', 'total_' ]:
                        assigned[f'{prefix}flux'] = assigned.apply(lambda row: flux_to_list(row, prefix, ''), axis=1)
                        assigned[f'{prefix}flux_err'] = assigned.apply(lambda row: flux_to_list(row, prefix, '_err'), axis=1)
                        
                    assigned['filter'] = assigned.apply(filter_to_list, axis=1)
                elif photometry is not None and photometry.bands is not None:
                    # Convert the flux columns into a single column of list objects
                    # Refer to canonical column names in __calculate_target_list_flux_bands

                    def flux_to_list(row, prefix='psf_flux_'):
                        fluxes = []
                        for b in photometry.bands.keys():
                            if row[f'filter_{b}'] is not None:
                                fluxes.append(row[f'{prefix}{b}'])
                        return fluxes
                    
                    def filter_to_list(row):
                        filters = []
                        for b in photometry.bands.keys():
                            if row[f'filter_{b}'] is not None:
                                filters.append(row[f'filter_{b}'])
                        return filters
                    
                    for prefix in [ 'psf_', 'fiber_', 'total_' ]:
                        assigned[f'{prefix}flux'] = assigned.apply(lambda row: flux_to_list(row, 'psf_flux_'), axis=1)
                        assigned[f'{prefix}flux_err'] = assigned.apply(lambda row: flux_to_list(row, 'psf_flux_err_'), axis=1)

                    assigned['filter'] = assigned.apply(filter_to_list, axis=1)
                else:
                    # No fluxes or magnitudes available (sky), just generate some empty lists
                    for prefix in [ 'psf_', 'fiber_', 'total_' ]:
                        assigned[f'{prefix}flux'] = assigned.apply(lambda row: [], axis=1)
                        assigned[f'{prefix}flux_err'] = assigned.apply(lambda row: [], axis=1)

                    assigned['filter'] = assigned.apply(lambda row: [], axis=1)

                assigned_target_lists[k] = assigned

        return assigned_target_lists
    
    def __merge_assigned_target_lists(self, assigned_target_lists):
        """
        Compile the list of the targets that will make it into the final design files.
        Some targets might be repeated in more than one target lists. These were taken care
        of in __append_source_target_lists by cross-matching each target list with the
        targets already added to the netflow target cache. Now we need to process the target
        lists in the same order and if a target is already in the assigned_targets_all list,
        we need to update the fluxes, magnitudes, etc. Otherwise, we append the target to the
        assigned_targets_all list.
        """

        # Take only a subset of columns, these will be used to generate the design file
        # NOTE: targetid is the index of each target list
        columns = [ '__target_idx', '__key',
                    'epoch', 'proposalid', 'tract', 'patch', 'catid', 'obcode',
                    'filter', 'psf_flux', 'psf_flux_err',
                    'fiber_flux', 'fiber_flux_err',
                    'total_flux', 'total_flux_err' ]
        
        assigned_targets_all = None
        for prefix in ['cal', 'sci']:
            for k, assigned_targets in assigned_target_lists.items():
                if self._config.targets[k].prefix == prefix:

                    if assigned_targets_all is None:
                        assigned_targets_all = assigned_targets[columns]
                    else:
                        # Some targets may be in multiple catalogs, so we need to join on __target_idx
                        # and update if there's a match or append if there's no match

                        # Find records in assigned_target_list that have matching __target_idx in assigned_targets_all
                        mask = assigned_targets['__target_idx'].isin(assigned_targets_all['__target_idx'])
                        
                        # Update the records in assigned_targets_all
                        # TODO: now it's done one-by-one but it could be sped up a little bit by using a join
                        at = assigned_targets.set_index('__target_idx')
                        at_all = assigned_targets_all.set_index('__target_idx')
                        for i in np.where(mask)[0]:
                            target_idx = assigned_targets.loc[i, '__target_idx']

                            # Get the reference to the list storing the filters in assigned_targets_all
                            filter_all = at_all.loc[target_idx, 'filter']

                            # Get the reference to lists store in assigned_targets_all
                            flux_all = {}
                            flux_err_all = {}
                            for flux_type in ['psf', 'fiber', 'total']:
                                if f'{flux_type}_flux' in at_all.columns:
                                    flux_all[flux_type] = at_all.loc[target_idx, f'{flux_type}_flux']
                                    flux_err_all[flux_type] = at_all.loc[target_idx, f'{flux_type}_flux_err']

                            # Loop over all filters in the target list and append the fluxes to the all targets'
                            # list if the filter doesn't exist yet
                            flux = {}
                            flux_err = {}
                            for flux_type in ['psf', 'fiber', 'total']:
                                if f'{flux_type}_flux' in at.columns:
                                    flux[flux_type] = at.loc[target_idx, f'{flux_type}_flux']
                                    flux_err[flux_type] = at.loc[target_idx, f'{flux_type}_flux_err']

                            for j, f in enumerate(at.loc[target_idx, 'filter']):
                                if f not in filter_all:
                                    filter_all.append(f)
                                    for flux_type in flux:
                                        # Append the fluxe and error
                                        flux_all[flux_type].append(flux[flux_type][j])
                                        flux_err_all[flux_type].append(flux_err[flux_type][j])

                        assigned_targets_all = at_all.reset_index()

                        # Append the rest of the rows
                        assigned_targets_all = pd.concat([ assigned_targets_all, assigned_targets[~mask][columns] ], ignore_index=True)

        # Sky positions are never cross-matched, so simply append them to the assigned targets list
        for k, assigned_targets in assigned_target_lists.items():
            if self._config.targets[k].prefix == 'sky':
                assigned_targets_all = pd.concat([ assigned_targets_all, assigned_targets[columns] ], ignore_index=True)
        
        # Find duplicate target_idx
        duplicates = assigned_targets_all[assigned_targets_all.duplicated('__target_idx', keep=False)]
        if duplicates.shape[0] > 0:
            raise NotImplementedError()

        return assigned_targets_all
    
    def __join_assignments(self, assignments, assigned_targets_all):
        """
        Joins in additional target information to the assignments list, such as
        fluxes, filter names, etc. that's not already in the assignments list.

        Parameters
        ----------

        assignments : DataFrame
            A list of fibers for each visit
        assigned_targets_all : DataFrame
            A merged list of targets that were assigned to any of the fibers
        """
        
        # Only include columns that are not already in `assignments`
        columns = assigned_targets_all.columns.difference(assignments.columns)

        # Join the assignments with the assigned targets to get the final object lists
        assignments_all = pd.merge(
            pd_to_nullable(assignments.set_index('target_idx')), 
            pd_to_nullable(assigned_targets_all[columns].set_index('__target_idx')),
            how='left', left_index=True, right_index=True)
        
        assignments_all.reset_index(inplace=True, names='__target_idx')

        return assignments_all
    
    def __substitute_nans(self, df):
        """
        Substitute NaN values with zeros or other default values.
        """

        special = {
            'parallax': 1e-7,
        }

        for c in df.columns:
            if c in special:
                value = special[c]
            else:
                value = 0.0

            pd_fillna(df, c, value)
            
    def __create_designs(self, netflow : Netflow, assignments_all):
        designs = []

        for visit in netflow.visits:
            design_name = f'{self._config.field.name} P{visit.pointing_idx:03d}/{len(netflow.pointings):03d} V{visit.visit_idx:02d}/{len(netflow.visits):02d}'
            d = Design.create_pfsDesign_visit(visit, assignments_all,
                                              design_name=design_name,
                                              arms=self._config.field.arms)
            designs.append(d)

            logger.info(f'Generated design {d.pfsDesignId:016x} for visit {visit.visit_idx} with name `{d.designName}`.')

        unique_id = np.unique([ d.pfsDesignId for d in designs])
        logger.info(f'Generated {len(designs)} designs from which {unique_id.size} are unique.')

        return designs
    
    def __verify_designs(self, designs):
        for design in designs:

            if not np.all(np.bincount(design.fiberId) <= 1):
                logger.error(f'Duplicate fiber assignments in design {design.pfsDesignId:016x}.')

            if not np.all(np.isnan(design.dec[design.targetType == TargetType.UNASSIGNED])):
                logger.error(f'Unassigned targets with coordinates found in design {design.pfsDesignId:016x}.')

            # Further verifications - probably better done in the notebook or in a separate script
            # - make sure unassigned cobras are in home position
            # - make sure all assigned fibers have flux given

            # design.pfiNominal[design.targetType == TargetType.UNASSIGNED]
            # 'na': TargetType.UNASSIGNED,
            # 'sky': TargetType.SKY,
            # 'cal': TargetType.FLUXSTD,
            # 'sci': TargetType.SCIENCE,
            # 'eng': TargetType.ENGINEERING,

    
    def __save_designs(self, designs):
        for d in designs:
            d.write(dirName=self.__outdir)

            # Make sure that the file is created
            fn = os.path.join(self.__outdir, PfsDesign.fileNameFormat % (d.pfsDesignId))
            if not os.path.isfile(fn):
                raise RuntimeError('Design file was not be written.')

            logger.info(f'Saved design {d.pfsDesignId:016x} to `{d.filename}`.')

    def _get_design_list_path(self):
        return os.path.join(self.__outdir, f'{self._config.field.key}_designs.feather')

if __name__ == '__main__':
    script = NetflowScript()
    script.execute()