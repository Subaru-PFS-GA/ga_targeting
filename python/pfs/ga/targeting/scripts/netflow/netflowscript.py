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

        self.__resume = False
        self.__stage = None
        self.__nvisits = None
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
        self.add_arg('--in', type=str, required=True, help='Path to the input directory with target lists.')
        self.add_arg('--out', type=str, required=True, help='Path to the output directory.')
        self.add_arg('--resume', action='store_true', help='Resume from a previous run.')
        self.add_arg('--stage', type=int, nargs='*', help='Stage of observations.')
        self.add_arg('--nvisits', type=int, help='Number of visits for each pointing.')
        self.add_arg('--exp-time', type=float, help='Exposure time per visit, in seconds.')
        self.add_arg('--obs-time', type=str, help='Observation time in ISO format in UTC.')
        self.add_arg('--time-limit', type=int, help='Time limit for the optimization in seconds.')
        
        self.add_arg('--xmatch-rad', type=float, help='Cross-match radius in arcseconds.')

        self.add_arg('--skip-notebooks', action='store_true', help='Skip execution of evaluation notebooks.')

    def _init_from_args_pre_logging(self, args):
        super()._init_from_args_pre_logging(args)

        self._indir = self.get_arg('in', args, self._indir)
        self._outdir = self.get_arg('out', args, self._outdir)
        self.__resume = self.get_arg('resume', args, self.__resume)

    def _init_from_args(self, args):
        super()._init_from_args(args)

        # If a field is specified, load its default configuration      
        self._create_field_from_args(args)
        self._create_config_from_field()

        # Load the configuration template files and merge with the default config
        self._load_config_files(args)

        self.__stage = self.get_arg('stage', args, self.__stage)
        self.__nvisits = self.get_arg('nvisits', args, self.__nvisits)
        self.__exp_time = self.get_arg('exp_time', args, self.__exp_time)
        self.__obs_time = self.get_arg('obs_time', args, self.__obs_time)
        self.__time_limit = self.get_arg('time_limit', args, self.__time_limit)
        self.__xmatch_rad = self.get_arg('xmatch_rad', args, self.__xmatch_rad)
        self.__skip_notebooks = self.get_arg('skip_notebooks', args, self.__skip_notebooks)

        # Override the configuration with the command-line arguments
        if self.__nvisits is not None:
            self._config.field.nvisits = self.__nvisits

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
        if not self.__resume and os.path.isdir(self._outdir):
            raise Exception(f'Output directory already exists: {self._outdir}')
        elif self.__resume and not os.path.isdir(self._outdir):
            raise Exception(f'Output directory does not exist: {self._outdir}, cannot resume.')
        elif not self.__resume:
            os.makedirs(self._outdir)

        # Update log file path to point to the output directory
        self.log_file = os.path.join(self._outdir, os.path.basename(self.log_file))

    def run(self):
        # Load instrument calibration data
        instrument = self._create_instrument()

        # Generate the list of pointings
        pointings = self._generate_pointings(stage=self.__stage)

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
        target_lists = self.load_preprocessed_target_lists()

        # Load the netflow target list from the input directory and run
        # the validation step again
        fn = os.path.join(self._indir, 'netflow_targets.feather')
        netflow.load_targets(fn)
        netflow.validate_targets()

        # Run the netflow optimization                
        self.__run_netflow(netflow)

        # TODO update netflow.__targets with done_visits
        #      then save netflow.__targets
        #      check if req_visits or done_visits are ever copied from the
        #      target lists or it's just exp_time that is used

        # Save the targets lists to the output directory
        # This is necessary only if we want to update anything after running
        # netflow but also because this way we can always run the next stage
        # from a single input directory
        self._save_preprocessed_target_lists(target_lists, self._outdir)

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
            for notebook in ['targets', 'assignments', 'calibration', 'cobra_groups', 'design']:
                logger.info(f'Executing evaluation notebook `{notebook}`...')
                notebook_path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), f'scripts/netflow/notebooks/{notebook}.ipynb')
                parameters = {
                    'DEBUG': False,
                    'CONFIG_FILE': self._get_output_config_path(),
                    'OUTPUT_PATH': self._outdir,
                }
                self._execute_notebook(notebook_path, parameters, self._outdir)
    
    def load_preprocessed_target_lists(self):
        target_lists = {}
        for k, target_list_config in self._config.targets.items():
            # Load the preprocessed target lists from the input directory
            # These will be later update with the number of scheduled visits
            # and saved to the output directory
            path = self._get_preprocessed_target_list_path(k, self._indir)
            target_lists[k] = self._load_preprocessed_target_list(k, path)

        return target_lists
    
    def __run_netflow(self, netflow):
        logger.info('Initializing netflow.')
        netflow.init(resume=self.__resume, save=True)

        logger.info('Building the netflow model.')
        netflow.build(resume=self.__resume, save=True)

        logger.info('Solving the netflow model.')
        netflow.solve()

        logger.info('Saving the netflow solution.')
        netflow.save_solution()

        logger.info('Extracting assignments and updating the target list')
        netflow.extract_assignments()
        netflow.update_done_visits()

        logger.info('Saving target list with updated done visit counts.')
        netflow.save_targets()

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
        return os.path.join(self._outdir, f'{self._config.field.key}_assignments.feather')
    
    def __save_assignments(self, assignments):
        path = self.__get_assignments_path()
        logger.info(f'Saving fiber assignments to `{path}`.')
        assignments.to_feather(path)

    def __load_assignments(self):
        path = self.__get_assignments_path()
        logger.info(f'Loading fiber assignments from `{path}`.')
        return pd.read_feather(path)
    
    def __get_summary_path(self):
        return os.path.join(self._outdir, f'{self._config.field.key}_summary.feather')
    
    def __save_summary(self, summary):
        path = self.__get_summary_path()
        logger.info(f'Saving target assignment summary to `{path}`.')
        summary.to_feather(path)
    
    def _get_assignments_all_path(self):
        return os.path.join(self._outdir, f'{self._config.field.key}_assignments_all.feather')

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
            d.write(dirName=self._outdir)

            # Make sure that the file is created
            fn = os.path.join(self._outdir, PfsDesign.fileNameFormat % (d.pfsDesignId))
            if not os.path.isfile(fn):
                raise RuntimeError('Design file was not be written.')

            logger.info(f'Saved design {d.pfsDesignId:016x} to `{d.filename}`.')

    def _get_design_list_path(self):
        return os.path.join(self._outdir, f'{self._config.field.key}_designs.feather')

if __name__ == '__main__':
    script = NetflowScript()
    script.execute()