import os
import commentjson as json
from datetime import datetime, timedelta, tzinfo
import pandas as pd
from pandas import Float32Dtype, Float64Dtype, Int32Dtype, Int64Dtype
import astropy.units as u

from ..config import NetflowConfig
from ..targets.dsph import GALAXIES
from ..instrument import SubaruPFI
from ..io import ObservationSerializer
from ..netflow import Netflow, Design
from ..core import Pointing
from ..util.astro import *
from ..util.pandas import *
from .script import Script

from ..setup_logger import logger

class NetflowScript(Script):
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
        * epoch: epoch of the coordinates, default is J2000.0. It has no effect if
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
        self.__nvisits = None
        self.__exp_time = None
        self.__obs_time = None
        self.__time_limit = None
        self.__nan_values = True
        
        self.__field = None
        self.__config = None

    def __get_config(self):
        return self.__config
    
    config = property(__get_config)

    def _add_args(self):
        super()._add_args()

        self.add_arg('--dsph', type=str, help='Name of a predefined target.')
        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')
        self.add_arg('--resume', action='store_true', help='Resume from a previous run.')
        self.add_arg('--nvisits', type=int, help='Number of visits for each pointing.')
        self.add_arg('--exp-time', type=float, help='Exposure time per visit, in seconds.')
        self.add_arg('--obs-time', type=str, help='Observation time in ISO format in UTC.')
        self.add_arg('--time-limit', type=int, help='Time limit for the optimization in seconds.')

        self.add_arg('--nan-values', action='store_true', dest='nan_values', help='Use NaN values in the output files.')
        self.add_arg('--no-nan-values', action='store_false', dest='nan_values', help='Do not use NaN values in the output files.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self.__outdir = self.get_arg('out', args, self.__outdir)
        self.__resume = self.get_arg('resume', args, self.__resume)

    def _init_from_args(self, args):
        super()._init_from_args(args)

        if self.is_arg('dsph', args):
            self.__field = GALAXIES[self.get_arg('dsph', args)]

        # If a field is specified, load its default configuration      
        if self.__field is not None:
            self.__config = self.__field.get_netflow_config()
        else:
            self.__config = NetflowConfig.default()

        # Load the configuration template files and merge with the default config
        config_files = self.get_arg('config', args)
        self.__config.load(config_files, ignore_collisions=True)
        logger.info(f'Loaded {len(config_files)} config files from {config_files}.')
        logger.info(f'Found {len(self.__config.targets)} target lists in the configuration.')

        self.__nvisits = self.get_arg('nvisits', args, self.__nvisits)
        self.__exp_time = self.get_arg('exp_time', args, self.__exp_time)
        self.__obs_time = self.get_arg('obs_time', args, self.__obs_time)
        self.__time_limit = self.get_arg('time_limit', args, self.__time_limit)
        self.__nan_values = self.get_arg('nan_values', args, self.__nan_values)

        # Override the configuration with the command-line arguments
        if self.__time_limit is not None:
            self.__config.gurobi_options.timelimit = self.__time_limit

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

    def _dump_settings(self):
        super()._dump_settings()

        # Save the active configuration to the output directory
        command = self.get_command_name()
        path = os.path.join(self.__outdir, f'{command}_{self.timestamp}.config')
        self.__config.save(path)

        logger.debug(f'Configuration saved to `{os.path.abspath(path)}`.')

    def run(self):

        # Load instrument calibration data
        instrument = self.__create_instrument()

        # Generate the list of pointings
        pointings = self.__generate_pointings()

        # Create the netflow object
        netflow = Netflow(
            f'{self.__config.field.name}',
            instrument,
            pointings,
            workdir = self.__outdir,
            netflow_options = self.__config.netflow_options,
            solver_options = self.__config.gurobi_options,
            debug_options = self.__config.debug_options)

        # Load the original target list files or, if resuming the processing,
        # load the preprocessed target lists from the output directory.
        target_lists = self.load_source_target_lists()

        # Look for duplicates, etc.
        self.__validate_source_target_lists(target_lists)

        # Append the source target lists to the netflow object
        self.__append_source_target_lists(netflow, target_lists)
            
        # TODO: plot the sky with PFI footprints
            
        self.__run_netflow(netflow)

        # We do not implement resume logic beyond this point because it's not worth it time-wise.

        # TODO: verify collisions and unassign cobras with lower priority targets
        self.__unassign_colliding_cobras(netflow)
        
        # Extract the assignments from the netflow solution
        assignments = self.__extract_assignments(netflow)
        self.__save_assignments(assignments)

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

        assignments_all = self.__join_assignments(assignments, assigned_targets_all)
        self.__save_assignments_all(assignments_all)

        # What's left is to write out the design files and the assigned targets
        # in the web uploader format. For this, substitute all NA and NaN values with
        # zeros
        if not self.__nan_values:
            self.__substitute_nans(assignments_all)

        # Generate the designs and save them to the output directory
        designs = self.__create_designs(netflow, assignments_all)

        # # Verify each design
        # # TODO: move this after any other step and throw exception on
        # #       verification failure, such az endpoint collisions or 
        # #       elbow collisions. Only report trajectory collisions.
        # logger.info('Verifying solution.') 
        # netflow.verify()

        self.__save_designs(designs)

        # Save the assigned targets in the official web uploader format
        for visit in netflow.visits:
            assignments_web = self.__merge_assignments_web(assignments, target_lists, visit.visit_idx)
            self.__save_assignments_web(assignments_web, visit.visit_idx, format='feather')
            self.__save_assignments_web(assignments_web, visit.visit_idx, format='csv')

    def __create_instrument(self):
        return SubaruPFI(instrument_options=self.__config.instrument_options)
    
    def load_source_target_lists(self):
        target_lists = {}
        for k, target_list_config in self.__config.targets.items():
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
    def load_target_list(key, target_list_config):
        if target_list_config.reader is not None:
            # TODO: implement custom reader function
            raise NotImplementedError()
        
        fn = os.path.expandvars(target_list_config.path)
        logger.info(f'Reading target list from `{fn}`.')

        reader = ObservationSerializer(
            catalog_name=key,
            columns = target_list_config.columns,
            column_map = target_list_config.column_map,
            data_types = target_list_config.data_types,
            value_map = target_list_config.value_map,
            index = target_list_config.index,
            filters = target_list_config.filters,
            bands = target_list_config.bands,
            limits = target_list_config.limits,
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

        if self.__config.targets[key].filters is not None:
            return self.__calculate_target_list_flux_filters(key, target_list)
        elif self.__config.targets[key].bands is not None:
            return self.__calculate_target_list_flux_bands(key, target_list)
        else:
            # No fluxes or magnitudes available (sky)
            pass

    def __calculate_target_list_flux_filters(self, key, target_list):
        """
        Calculate the missing flux values from other columns that are available
        """

        target_list.calculate_flux_filters(self.__config.targets[key].filters)

    def __calculate_target_list_flux_bands(self, key, target_list):
        """
        Calculate the missing flux values from other columns that are available
        """

        target_list.calculate_flux_bands(self.__config.targets[key].bands)

    def __append_target_list_extra_columns(self, key, target_list):
        """
        Join the assignments with the target lists to append the fluxes that will
        be included in the design files.
        """

        logger.info(f'Appending extra columns to target list `{key}`.')

        # Substitute these fields into column patterns
        field_info = dict(
            name=self.__config.field.key,               # Field short name, e.g. 'umi'
            key=key,                                    # Target list key, e.g. 'dsph' or 'fluxstd'
            obs_time=self.__config.field.obs_time,
            resolution=self.__config.field.resolution,
            prefix=self.__config.targets[key].prefix,
            catid=self.__config.targets[key].catid,
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
        if self.__config.targets[key].extra_columns is not None:
            for c, config in self.__config.targets[key].extra_columns.items():
                patterns = [config.constant, config.lambda_func, config.pattern]
                if np.sum([p is not None for p in patterns]) > 1:
                    raise Exception(f'Only one type of pattern is allowed for extra column `{c}`.')
                
            # The order of column generation is: constants, lambda, pattern
            for c, config in self.__config.targets[key].extra_columns.items():
                if config.constant is not None:
                    set_extra_column_constant(c, config.dtype, None, config.constant)
            
            for c, config in self.__config.targets[key].extra_columns.items():
                if config.lambda_func is not None:
                    set_extra_column_lambda(c, config.dtype, None, config.lambda_func, config.lambda_args)

            for c, config in self.__config.targets[key].extra_columns.items():
                if config.pattern is not None:
                    set_extra_column_pattern(c, config.dtype, None, config.pattern)

        # If constants are defined for these columns in the config, assign them
        set_extra_column_constant('epoch', 'string', 'J2000.0', self.__config.targets[key].epoch)
        set_extra_column_constant('catid', pd.Int32Dtype(), -1, self.__config.targets[key].catid)
        set_extra_column_constant('priority', pd.Int32Dtype(), None, self.__config.targets[key].priority)
        set_extra_column_constant('exp_time', pd.Int32Dtype(), None, self.__config.targets[key].exp_time)

        # Create a targetid column automatically if it's not present
        if 'targetid' not in target_list.columns:
            target_list.append_column('targetid', target_list.data.index, pd.Int64Dtype())
            logger.warning(f'Automatically generated column `targetid` in target list `{key}`.')

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
            should.extend(['pmra', 'pmdec', 'parallax', 'epoch'])
            
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

        logger.info(f'Successfully validated target list `{key}`.')

    def __validate_source_target_lists(self, target_lists):
        # TODO: make sure all targetids are unique across all target lists
        pass

    def __validate_fluxstd_targets(self, target_lists):
        # TODO: make sure cobra group constraints are satisfiable
        pass

    def __validate_sky_targets(self, target_lists):
        # TODO: make sure cobra group constraints are satisfiable
        pass

    def __append_source_target_lists(self, netflow, target_lists):
        for k, target_list in target_lists.items():
            if self.__config.targets[k].prefix == 'sci':
                netflow.append_science_targets(target_list)
            elif self.__config.targets[k].prefix == 'cal':
                netflow.append_fluxstd_targets(target_list)
                        
                # Make sure cobra location/instrument group constraints are satisfiable
                self.__validate_fluxstd_targets(target_lists)
            elif self.__config.targets[k].prefix == 'sky':
                netflow.append_sky_targets(target_list)

                # Make sure cobra location/instrument group constraints are satisfiable
                self.__validate_sky_targets(target_lists)
            else:
                raise NotImplementedError()
    
    def __get_preprocessed_target_list_path(self, key):
        return os.path.join(self.__outdir, f'{self.__config.field.key}_targets_{key}.feather')

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

    def __generate_pointings(self):

        # The list of pointings can be defined in the config file, which then
        # override the pointings defined for the field in the library source code.

        if self.__config.pointings is not None:
            pp = [ p.get_pointing() for p in self.__config.pointings ]
        elif self.__field is not None:
            # Load from the class
            pp = self.__field.get_pointings(SubaruPFI)
        else:
            raise NotImplementedError()
        
        # The number of visits can be overridden from the command-line argument
        nvisits = self.__nvisits if self.__nvisits is not None else self.__config.field.nvisits
        exp_time = self.__exp_time if self.__exp_time is not None else self.__config.field.exp_time

        if self.__obs_time is not None:
            obs_time = datetime.fromisoformat(self.__obs_time)
        else:
            obs_time = self.__config.field.obs_time
        

        # Create new pointing objects with obs_time, exp_time etc.
        pointings = []
        for p in pp:   
            pointings.append(Pointing(
                p.ra, p.dec,
                posang = p.posang,
                obs_time = obs_time,
                exp_time = exp_time,
                nvisits = nvisits))
            
        return pointings
    
    def __run_netflow(self, netflow):
        logger.info('Building netflow model.')
        netflow.build(resume=self.__resume, save=True)

        logger.info('Solving netflow model.')
        netflow.solve()

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
        return os.path.join(self.__outdir, f'{self.__config.field.key}_assignments.feather')
    
    def __save_assignments(self, assignments):
        path = self.__get_assignments_path()
        logger.info(f'Saving fiber assignments to `{path}`.')
        assignments.to_feather(path)

    def __load_assignments(self):
        path = self.__get_assignments_path()
        logger.info(f'Loading fiber assignments from `{path}`.')
        return pd.read_feather(path)
    
    def __get_summary_path(self):
        return os.path.join(self.__outdir, f'{self.__config.field.key}_summary.feather')
    
    def __save_summary(self, summary):
        path = self.__get_summary_path()
        logger.info(f'Saving target assignment summary to `{path}`.')
        summary.to_feather(path)
    
    def __get_assignments_all_path(self):
        return os.path.join(self.__outdir, f'{self.__config.field.key}_assignments_all.feather')
    
    def __save_assignments_all(self, assignments_all):
        # This dataframe contains columns with dtype `object`` that contain the list of
        # magnitudes, fluxes, filter names, etc.
        
        path = self.__get_assignments_all_path()
        logger.info(f'Saving assignments joined with the input catalogs to `{path}`.')
        assignments_all.to_feather(path)

    def __load_assignments_all(self):
        path = self.__get_assignments_all_path()
        logger.info(f'Loading assignments joined with the input catalogs from `{path}`.')
        return pd.read_feather(path)

    def __append_flux_filter_lists(self, assignments, target_lists):
        """
        Join the assignments with the target lists to append the fluxes that will
        be included in the design files.
        """

        # Sky fibers and engeering fibers have no targetid, so make a unique list of
        # actual targetids
        unique_targetid = pd.DataFrame({ 'targetid': assignments['targetid'].unique() })
        unique_targetid.set_index('targetid', inplace=True)

        assigned_target_lists = {}
        for k, target_list in target_lists.items():
            assigned = target_list.data.set_index('targetid').join(unique_targetid, how='inner')

            # Filter names are either defined in the config or can be extracted from the columns
            # If the config entry is a dictionary, we assume that the keys are the filter names.
            # If it is a list, we assume that the filter names are stored as values in a column.

            if self.__config.targets[k].filters is not None:
                # Convert the flux columns into a single column of list objects
                # Refer to canonical column names in __calculate_target_list_flux_filters

                def flux_to_list(row, prefix, postfix):    
                    return [ row[f'{b}_{prefix}flux{postfix}'] for b in self.__config.targets[k].filters.keys() ]
                
                def filter_to_list(row):
                    return list(self.__config.targets[k].filters.keys())
                
                for prefix in [ 'psf_', 'fiber_', 'total_' ]:
                    assigned[f'{prefix}flux'] = assigned.apply(lambda row: flux_to_list(row, prefix, ''), axis=1)
                    assigned[f'{prefix}flux_err'] = assigned.apply(lambda row: flux_to_list(row, prefix, '_err'), axis=1)
                    
                assigned['filter'] = assigned.apply(filter_to_list, axis=1)
            elif self.__config.targets[k].bands is not None:
                # Convert the flux columns into a single column of list objects
                # Refer to canonical column names in __calculate_target_list_flux_bands

                def flux_to_list(row, prefix='psf_flux_'):
                    fluxes = []
                    for b in self.__config.targets[k].bands.keys():
                        if row[f'filter_{b}'] is not None:
                            fluxes.append(row[f'{prefix}{b}'])
                    return fluxes
                
                def filter_to_list(row):
                    filters = []
                    for b in self.__config.targets[k].bands.keys():
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
        # Take only a subset of columns, these will be used to generate the design file
        # NOTE: targetid is the index of each target list
        columns = [ 'epoch', 'proposalid', 'tract', 'patch', 'catid', 'obcode',
                    'filter', 'psf_flux', 'psf_flux_err',
                    'fiber_flux', 'fiber_flux_err',
                    'total_flux', 'total_flux_err' ]
        
        assigned_targets_all = None
        for key, assigned_target_list in assigned_target_lists.items():
            if assigned_targets_all is None:
                assigned_targets_all = assigned_target_list[columns]
            else:
                assigned_targets_all = pd.concat([ assigned_targets_all, assigned_target_list[columns] ])

        return assigned_targets_all
    
    def __join_assignments(self, assignments, assigned_targets_all):
        # `assignments` is a list of fibers for every visit
        # `assigned_targets_all` is a list of targets that were assigned to any of the fibers

        # Only include columns that are not already in `assignments`
        columns = assigned_targets_all.columns.difference(assignments.columns)

        # Join the assignments with the assigned targets to get the final object lists
        assignments_all = pd.merge(
            pd_to_nullable(assignments.set_index('targetid')), 
            pd_to_nullable(assigned_targets_all[columns]),
            how='left', left_index=True, right_index=True)
        
        assignments_all.reset_index(inplace=True, names='targetid')

        # TODO: substitute NaN values

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

            if df.dtypes[c] == float or df.dtypes[c] == np.float64 or df.dtypes[c] == np.float32:
                df.loc[df[c].isna(), c] = value
            elif df.dtypes[c] == Float64Dtype() or df.dtypes[c] == Float32Dtype():
                df[c] = df[c].fillna(value)
                df.loc[df[c].isna(), c] = value
            elif df.dtypes[c] == 'string':
                df[c] = df[c].fillna('')
            elif df.dtypes[c] == object:
                pass
            else:
                pass
            
    def __create_designs(self, netflow, assignments_all):
        designs = []

        for visit in netflow.visits:
            d = Design.create_pfsDesign_visit(visit, assignments_all, arms=self.__config.field.arms)
            d.designName = f'ga_{self.__config.field.name}'
            designs.append(d)

            logger.info(f'Generated design {d.pfsDesignId:016x} for visit {visit.visit_idx} with name `{d.designName}`.')

        return designs
    
    def __save_designs(self, designs):
        for d in designs:
            d.write(dirName=self.__outdir)
            logger.info(f'Saved design {d.pfsDesignId:016x} to `{d.filename}`.')

    def __merge_assignments_web(self, assignments, target_lists, visit_idx):
        """
        Generate an assignment list that is compatible with the web uploader. This is for compatibility
        with the rest of the targeting algorithms.
        """

        assignment_columns = ['targetid', 'RA', 'Dec', 'pm', 'pmra', 'pmdec', 'parallax', 'exp_time', 'priority']
        target_list_columns = ['targetid', 'obcode']

        assignments_web = None
        all_flux_columns = set()

        for k, target_list in target_lists.items():
            
            # Collect all available flux columns from the target list
            column_map = {}         # column renames
            extra_columns = {}      # flux columns calculated from magnitudes and other extra columns
            
            if target_list.photometry is not None:
                for p, photometry in target_list.photometry.items():
                    for m, magnitude in photometry.magnitudes.items():
                        mag, mag_err, flux, flux_err, ext, magnitude_type = \
                            target_list.get_magnitude_column_names(magnitude)
                        
                        # Select the flux column and its error and format column names
                        # into the web uploader format. If the flux is not available,
                        # calculate it from the magnitude. The flux unit is nJy.
                        if flux is not None:
                            name = f'{m}_{p}'
                            extra_columns[name] = target_list.data[flux]
                            all_flux_columns.add(name)
                            if flux_err is not None:
                                name = f'{m}_{p}_error'
                                extra_columns[name] = target_list.data[flux_err]
                                all_flux_columns.add(name)
                        elif mag is not None:
                            mag = target_list.data[mag] if mag is not None and mag in target_list.data else None
                            mag_err = target_list.data[mag_err] if mag_err is not None and mag_err in target_list.data else None

                            flux, flux_err = ABmag_to_nJy(mag, mag_err)
                            if flux is not None:
                                name = f'{m}_{p}'
                                extra_columns[name] = flux
                                all_flux_columns.add(name)
                                if flux_err is not None:
                                    name = f'{m}_{p}_error'
                                    extra_columns[name] = flux_err
                                    all_flux_columns.add(name)

            mask = assignments['visit_idx'] == visit_idx
            a = pd_to_nullable(assignments[assignment_columns][mask])
            a = a.set_index('targetid')

            b = target_list.data[target_list_columns].assign(**extra_columns)
            b = pd_to_nullable(b).rename(columns=column_map).set_index('targetid')

            aw = pd.merge(a, b, how='inner', left_index=True, right_index=True)
            aw.reset_index(inplace=True, names='targetid')
            aw.rename(inplace=True,
                      columns={
                          'targetid': 'obj_id',
                          'obcode': 'ob_code',
                          'RA': 'ra',
                          'Dec': 'dec',
                          'exp_time': 'exptime',
                    })
            
            # TODO: make these command-line arguments
            aw['resolution'] = 'm'
            aw['reference_arm'] = 'r'

            # Substitute <N/A> values
            for c in ['pm', 'pmra', 'pmdec']:
                aw[c] = aw[c].fillna(0.0).astype('float64')
            for c in ['parallax']:
                aw[c] = aw[c].fillna(1e-7).astype('float64')
            
            # Add columns: tract, patch ??

            if assignments_web is None:
                assignments_web = aw
            else:
                assignments_web = pd.concat([assignments_web, aw])

        # Substitute remaining <N/A> values with NaN
        pd_null_to_nan(assignments_web, columns=all_flux_columns, in_place=True)

        # Reorder columns
        columns = ['obj_id', 'ob_code', 'ra', 'dec', 'pm', 'pmra', 'pmdec', 'parallax', 'exptime', 'priority', 'resolution', 'reference_arm']
        for c in assignments_web.columns:
            if c not in columns:
                columns.append(c)
        assignments_web = assignments_web.reindex(columns, axis=1)

        return assignments_web
    
    def __get_assignments_web_path(self, visit_idx, format='feather'):
        return os.path.join(self.__outdir, f'{self.__config.field.key}_assignments_web_{visit_idx:03d}.{format}')
    
    def __save_assignments_web(self, assignments_web, visit_idx, format='feather'):
        # This dataframe contains columns with dtype `object`` that contain the list of
        # magnitudes, fluxes, filter names, etc.
        
        path = self.__get_assignments_web_path(visit_idx, format=format)
        logger.info(f'Saving assignments for web upload to `{path}`.')
        if format == 'csv':
            assignments_web.to_csv(path, index=False)
        elif format == 'feather':
            assignments_web.to_feather(path)
        else:
            raise NotImplementedError()

if __name__ == '__main__':
    script = NetflowScript()
    script.execute()