import os
import commentjson as json
import pandas as pd
import astropy.units as u

from ..config import NetflowConfig
from ..instrument import SubaruPFI
from ..netflow import Netflow, Design
from ..netflow.util import *
from ..core import Pointing
from ..util.astro import *
from .script import Script

from ..setup_logger import logger

class Assign(Script):
    """
    Command-line script to execute netflow optimazion and generate
    the target list and design files.
    """

    def __init__(self):
        super().__init__()

        self.__outdir = None

        self.__config = None

    def _add_args(self):
        super()._add_args()

        self.add_arg('--config', type=str, required=True, nargs='+', help='Path to the configuration file.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')

    def _init_from_args(self, args):
        super()._init_from_args(args)

        # Read the configuration file, assuming it a python file
        self.__config = NetflowConfig()

        # Load the configuration template file
        config_files = self.get_arg('config', args)
        self.__config.load(config_files, ignore_collisions=True)

        self.__outdir = self.get_arg('out', args, self.__outdir)

    def prepare(self):
        super().prepare()

        # Create the output directory
        if os.path.isdir(self.__outdir):
            raise Exception(f'Output directory already exists: {self.__outdir}')
        else:
            os.makedirs(self.__outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self.__outdir, os.path.basename(self.log_file))

        # Save the active configuration to the output directory
        command = self.get_command_name()
        self.__config.save(os.path.join(self.__outdir, f'{command}_{self.timestamp}.config.json'))

    def run(self):
        logger.info(f'Found {len(self.__config.targets)} target lists in the configuration.')

        target_lists = {}
        for k, target_list_config in self.__config.targets.items():
            target_list = self.__load_target_list(target_list_config)
            self.__validate_target_list(k, target_list_config, target_list)
            self.__calculate_target_list_flux(k, target_list)
            self.__append_target_list_extra_columns(k, target_list)
            target_lists[k] = target_list

        self.__validate_all_target_lists(target_lists)

        # Save the preprocessed target lists to the output direcotory
        for k, target_list_config in self.__config.targets.items():
            self.__save_target_list(k, target_lists[k])

        # TODO: plot the sky with PFI footprints
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
        
        # TODO: add functions to append targets as DataFrame,
        #       without using Catalog or Observation classes
        
        # Append the targets
        for k, target_list in target_lists.items():
            if self.__config.targets[k].prefix == 'sci':
                netflow.append_science_targets(target_list)
            elif self.__config.targets[k].prefix == 'cal':
                netflow.append_fluxstd_targets(target_list)
            elif self.__config.targets[k].prefix == 'sky':
                netflow.append_sky_targets(target_list)
            else:
                raise NotImplementedError()
            
        # TODO: plot the sky with PFI footprints

        for prefix in ['sci', 'cal', 'sky']:
            self.__save_netflow_targets(netflow, prefix)
            
        self.__run_netflow(netflow)
        
        # Extract the assignments from the netflow solution
        assignments = self.__extract_assignments(netflow)
        self.__save_assignments(assignments)

        # TODO: plot assignment statistics

        # Join the targets lists with the assignments to append the fluxes and
        # filter names as columns of lists, as  required for the design files
        assigned_target_lists = self.__append_flux_filter_lists(assignments, target_lists)

        # Concatenate the assigned target lists into a single DataFrame with only
        # the columns that are needed for the design files
        assigned_targets_all = self.__merge_assigned_target_lists(assigned_target_lists)

        # TODO: this is the point to postfix obcode with the visit id
        assignments_all = self.__join_assignments(assignments, assigned_targets_all)

        # Generate the designs and save them to the output directory
        designs = self.__create_designs(netflow, assignments_all)
        self.__save_designs(designs)

    def __create_instrument(self):
        return SubaruPFI(instrument_options=self.__config.instrument_options)

    def __load_target_list(self, target_list_config):
        if target_list_config.reader is not None:
            # TODO: implement custom reader function
            raise NotImplementedError()
        
        logger.info(f'Reading target list from `{target_list_config.path}`.')
        
        dir, filename = os.path.split(target_list_config.path)
        basename, ext = os.path.splitext(filename)

        if ext == '.feather':
            df = pd.read_feather(target_list_config.path, **target_list_config.reader_args)
        elif ext == '.csv':
            df = pd.read_csv(target_list_config.path, **target_list_config.reader_args)

        # If the column mapping is not None, rename columns accordingly
        if target_list_config.column_map is not None:
            df = df.rename(columns=target_list_config.column_map)

        logger.info(f'Read {len(df)} targets from target list.')

        return df
    
    def __save_target_list(self, key, target_list):
        fn = os.path.join(self.__outdir, f'{self.__config.field.name}_targets_{key}.feather')
        logger.info(f'Saving target list to `{fn}`.')
        target_list.to_feather(fn)
    
    def __calculate_target_list_flux(self, key, target_list):
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

        for filter, columns in self.__config.targets[key].filters.items():
            any_flux = None         # First available flux value for a filter
            any_flux_err = None     # First available flux error for a filter
            any_flux_key = None     # Key of the first available flux column
            for prefix in ['', 'psf_', 'fiber_', 'total_']:
                flux_key = prefix + 'flux'
                flux_err_key = prefix + 'flux_err'
                mag_key = prefix + 'mag'
                mag_err_key = prefix + 'mag_err'

                flux_col_canonical = f'{filter}_{prefix}flux'
                flux_err_col_canonical = f'{filter}_{prefix}flux_err'

                flux_col = columns[flux_key] if flux_key in columns else None
                flux_err_col = columns[flux_err_key] if flux_err_key in columns else None
                mag_col = columns[mag_key] if mag_key in columns else None
                mag_err_col = columns[mag_err_key] if mag_err_key in columns else None

                # Calculate flux from the magnitude
                flux = None
                flux_err = None
                if flux_col is None:
                    logger.warning(f'Missing flux column `{flux_key}` in target list {key}.')

                    if mag_col is not None:
                        # The magnitude is available to calculate the flux from
                        logger.info(f'Calculating `{flux_key}` from `{mag_key}` in filter {filter}.')

                        mag = target_list[mag_col]
                        mag_err = target_list[mag_err_col] if mag_err_col is not None else None

                        flux = (np.array(mag) * u.ABmag).to_value(u.nJy)
                        if mag_err is not None:
                            flux_err = 0.4 * np.log(10) * flux * np.array(mag_err)
                    elif any_flux is not None:
                        # No magnitude is available, copy flux with another prefix
                        logger.info(f'Copying `{flux_key}` from `{any_flux_key}` in filter {filter}.')
                        flux = any_flux
                        flux_err = any_flux_err
                else:
                    flux = target_list[flux_col]
                    if flux_err_col is not None:
                        flux_err = target_list[flux_err_col]

                # Save the flux to the target list with canonical column name
                target_list[flux_col_canonical] = flux if flux is not None else np.nan
                target_list[flux_err_col_canonical] = flux_err if flux_err is not None else np.nan

                if any_flux is None:
                    any_flux = target_list[flux_col_canonical]
                    any_flux_err = target_list[flux_err_col_canonical]
                    any_flux_key = flux_key

    def __calculate_target_list_flux_bands(self, key, target_list):

        for band, columns in self.__config.targets[key].bands.items():
            any_flux = None
            any_flux_err = None
            any_flux_key = None
            for prefix in ['', 'psf_', 'fiber_', 'total_']:
                # Column keys
                flux_key = prefix + 'flux'
                flux_err_key = prefix + 'flux_err'
                mag_key = prefix + 'mag'
                mag_err_key = prefix + 'mag_err'

                # Column names
                flux_col = columns[flux_key] if flux_key in columns else None
                flux_err_col = columns[flux_err_key] if flux_err_key in columns else None
                mag_col = columns[mag_key] if mag_key in columns else None
                mag_err_col = columns[mag_err_key] if mag_err_key in columns else None

                flux_col_canonical = f'{prefix}flux_{band}'
                flux_err_col_canonical = f'{prefix}flux_err_{band}'

                # Calculate flux from the magnitude
                flux = None
                flux_err = None
                if flux_col is None:
                    logger.warning(f'Missing flux column `{flux_key}` in target list {key}.')

                    if mag_col is not None:
                        # The magnitude is available to calculate the flux from
                        logger.info(f'Calculating `{flux_key}` from `{mag_key}` in filter {filter}.')

                        mag = target_list[mag_col]
                        mag_err = target_list[mag_err_col] if mag_err_col is not None else None

                        flux = (np.array(mag) * u.ABmag).to_value(u.nJy)
                        if mag_err is not None:
                            flux_err = 0.4 * np.log(10) * flux * np.array(mag_err)
                    elif any_flux is not None:
                        # No magnitude is available, copy flux with another prefix
                        logger.info(f'Copying `{flux_key}` from `{any_flux_key}` in filter {filter}.')
                        flux = any_flux
                        flux_err = any_flux_err
                else:
                    flux = target_list[flux_col]
                    if flux_err_col is not None:
                        flux_err = target_list[flux_err_col]

                # Save the flux to the target list with canonical column name
                target_list[flux_col_canonical] = flux if flux is not None else np.nan
                target_list[flux_err_col_canonical] = flux_err if flux_err is not None else np.nan

                if any_flux is None:
                    any_flux = target_list[flux_col_canonical]
                    any_flux_err = target_list[flux_err_col_canonical]
                    any_flux_key = flux_key

    def __append_target_list_extra_columns(self, key, target_list):
        """
        Join the assignments with the target lists to append the fluxes that will
        be included in the design files.
        """

        # If constants are defined for these columns in the config, assign them
        if self.__config.targets[key].epoch is not None:
            pd_append_column(target_list, 'epoch', self.__config.targets[key].epoch, 'string')
        elif 'epoch' not in target_list.columns:
            pd_append_column(target_list, 'epoch', 'J2000.0', 'string')
        else:
            target_list['epoch'] = target_list['epoch'].astype('string')

        if self.__config.targets[key].catid is not None:    
            pd_append_column(target_list, 'catid', self.__config.targets[key].catid, pd.Int32Dtype())
        elif 'catid' not in target_list.columns:
            pd_append_column(target_list, 'catid', -1, pd.Int32Dtype())
        else:
            target_list['catid'] = target_list['catid'].astype(pd.Int32Dtype())

        if self.__config.targets[key].proposalid is not None:
            pd_append_column(target_list, 'proposalid', self.__config.targets[key].proposalid, 'string') 
        elif 'proposalid' not in target_list.columns:
            pd_append_column(target_list, 'proposalid', '', 'string')
        else:
            target_list['proposalid'] = target_list['proposalid'].astype('string')

        # TODO: these must be calculated from the coordinates
        if 'tract' not in target_list.columns:
            pd_append_column(target_list, 'tract', 0, pd.Int32Dtype())
        else:
            target_list['tract'] = target_list['tract'].astype(pd.Int32Dtype())
        
        if 'patch' not in target_list.columns:
            pd_append_column(target_list, 'patch', '0,0', 'string')
        else:
            target_list['patch'] = target_list['patch'].astype('string')

        # TODO: how do we generate the obcodes?
        if 'obcode' not in target_list.columns:
            target_list['obcode'] = target_list['targetid'].astype('string')
        else:
            target_list['obcode'] = target_list['obcode'].astype('string')
    
    def __validate_target_list(self, key, target_list_config, target_list):
        # Verify that columns are present

        logger.info(f'Validating target list {key}.')

        must = ['targetid', 'RA', 'Dec']
        if target_list_config.prefix == 'sci':
            must.extend(['priority', 'exp_time'])

        for col in must:
            if col not in target_list.columns:
                raise Exception(f'Missing column "{col}" in target list {key}.')

        should = []
        if target_list_config.prefix in ['sci', 'cal']:
            should.extend(['pmra', 'pmdec', 'parallax', 'epoch'])
            
        for col in should:
            if col not in target_list.columns:
                logger.warning(f'Missing column "{col}" in target list {key}.')

        # Verify that targetid is unique and print duplicates
        duplicates = target_list[target_list.duplicated('targetid', keep=False)]
        if not duplicates.empty:
            logger.warning(f'Duplicate targetids found in target list {key}.')
            logger.warning(list(duplicates['targetid']))

        # TODO: verify that all columns are present that are in the filter list
        #       epoch, catid, proposalid, etc.

        logger.info(f'Successfully validated target list {key}.')

    def __validate_all_target_lists(self, target_lists):
        # TODO: make sure all targetids are unique across all target lists
        pass

    def __generate_pointings(self):
        if self.__config.field.pointings is not None:
            pp = [ p.get_pointing() for p in self.__config.field.pointings ]
        elif self.__config.field.definition is not None:
            # Load from the class
            pp = self.__config.field.definition.get_pointings(SubaruPFI)
        else:
            raise NotImplementedError()

        # Create new pointing objects with obs_time, exp_time etc.
        pointings = []
        for p in pp:   
            pointings.append(Pointing(
                p.ra, p.dec,
                posang=p.posang,
                obs_time=self.__config.field.obs_time,
                exp_time=self.__config.field.exp_time,
                nvisits=self.__config.field.nvisits))
            
        return pointings
    
    def __save_netflow_targets(self, netflow, prefix):
        fn = os.path.join(self.__outdir, f'netflow_targets_{prefix}.feather')
        logger.info(f'Saving targets with prefix `{prefix}` to `{fn}`.')
        mask = netflow.targets['prefix'] == prefix
        netflow.targets[mask].to_feather(fn)
    
    def __run_netflow(self, netflow):
        logger.info('Building netflow model.')
        netflow.build()
        
        logger.info('Saving netflow model.')
        netflow.save_problem()

        logger.info('Solving netflow model.')
        netflow.solve()

        logger.info('Saving netflow solution.')
        netflow.save_solution()
        
        logger.info('Extracting fiber assignments.')

    def __extract_assignments(self, netflow):
        """
        Extract the assignments from the netflow problem as a DataFrame with all
        necessary information to generate the design files.
        """
        
        assignments = netflow.get_target_assignments(
            include_target_columns=True,
            include_unassigned_fibers=True,
            include_engineering_fibers=True)

        return assignments
    
    def __save_assignments(self, assignments):
        fn = os.path.join(self.__outdir, f'{self.__config.field.name}_assignments.feather')
        logger.info(f'Saving fiber assignments to `{fn}`.')
        assignments.to_feather(fn)

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
            assigned = target_list.set_index('targetid').join(unique_targetid, how='inner')

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

        # Substitute NaN values


        return assignments_all
            
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

if __name__ == '__main__':
    assign = Assign()
    assign.execute()