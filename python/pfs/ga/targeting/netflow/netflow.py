import os
import logging
from typing import Callable
from collections import defaultdict
import numpy as np
import pandas as pd
from types import SimpleNamespace  

from ics.cobraOps.Bench import Bench

from ..util.args import *
from ..projection import Pointing
from ..instrument import SubaruWFC, SubaruPFI
from .calibtarget import CalibTarget
from .sciencetarget import ScienceTarget
from .gurobiproblem import GurobiProblem

class Netflow():
    # Wraps the netflow library with a more object oriented interface
    
    def __init__(self, 
                 name,
                 visits=None,
                 workdir=None,
                 filename_prefix=None,
                 solver_options=None,
                 netflow_options=None,
                 debug_options=None):

        # Configurable options

        self.__name = name                          # Problem name
        self.__visits = visits                      # List of pointings
        
        self.__workdir = workdir if workdir is not None else os.getcwd()
        self.__filename_prefix = filename_prefix if filename_prefix is not None else ''
        self.__solver_options = solver_options
        self.__netflow_options = netflow_options
        self.__debug_options = debug_options

        # Internally used variables

        self.__name_counter = 0

        self.__telescope_type = SubaruWFC
        self.__positioner_type = SubaruPFI
        self.__positioner = None

        self.__problem_type = GurobiProblem
        self.__problem = None                       # ILP problem, already wrapped

        self.__variables = None
        self.__constraints = None

        self.__target_classes = None
        self.__forbidden_pairs = None
        self.__black_dots = None
        self.__cobra_groups = None
        self.__time_budgets = None

        self.__visit_exp_time = None                # Exposure time of a single visit in integer seconds
        self.__targets = None                       # DataFrame of targets
        self.__target_fp_pos = None                 # Target focal plane positions for each pointing
        
        self.__assignments = None

    def __get_name(self):
        return self.__name
    
    def __set_name(self, value):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_targets(self):
        return self.__targets
    
    targets = property(__get_targets)

    def __get_target_classes(self):
        return self.__target_classes
    
    target_classes = property(__get_target_classes)

    def __get_assignments(self):
        return self.__assignments
    
    assignments = property(__get_assignments)

    def __get_netflow_option(self, key, default=None):
        if self.__netflow_options is not None and key in self.__netflow_options:
            return self.__netflow_options[key]
        else:
            return default
        
    def __get_debug_option(self, key, default=None):
        if self.__debug_options is not None and key in self.__debug_options:
            return self.__debug_options[key]
        else:
            return default

    def __get_prefixed_filename(self, filename=None, default_filename=None):
        # TODO: add support for .zip, .gz, .bz2, and .7z

        return os.path.join(
            self.__workdir,
            self.__filename_prefix,
            filename if filename is not None else default_filename
        )
    
    def __append_extension(self, filename, ext):
        # TODO: add support for .zip, .gz, .bz2, and .7z

        if os.path.splitext(filename)[1] != ext:
            filename += ext
        return filename
    
    def __get_target_filename(self, filename=None):
        # TODO: add support for .zip, .gz, .bz2, and .7z
        fn = self.__get_prefixed_filename(filename, 'netflow_targets.csv')
        fn = self.__append_extension(fn, '.csv')
        return fn
    
    def __get_problem_filename(self, filename=None):
        # TODO: add support for .zip, .gz, .bz2, and .7z
        fn = self.__get_prefixed_filename(filename, 'netflow_problem.mps')
        fn = self.__append_extension(fn, '.mps')
        return fn
    
    def __get_solution_filename(self, filename=None):
        # TODO: add support for .zip, .gz, .bz2, and .7z
        fn = self.__get_prefixed_filename(filename, 'netflow_solution.sol')
        fn = self.__append_extension(fn, '.sol')
        return fn
    
    def save_targets(self, filename=None):
        """Save the target list to a file"""

        raise NotImplementedError()

    def load_targets(self, filename=None):
        """Load the target list from a file"""

        raise NotImplementedError()

    def load_problem(self, filename=None, options=None):
        """Read the LP problem from a file"""

        fn = self.__get_problem_filename(filename)
        self.__problem = self.__problem_type(extraOptions=options)
        self.__problem.read_problem(fn)

    def save_problem(self, filename=None):
        """Save the LP problem to a file"""
        
        fn = self.__get_problem_filename(filename)
        self.__problem.write_problem(fn)
        
    def load_solution(self, filename=None):
        """Load an LP solution from a file."""

        fn = self.__get_solution_filename(filename)
        self.__problem.read_solution(fn)

    def save_solution(self, filename=None):
        """Save the LP solution to a file."""

        fn = self.__get_solution_filename(filename)
        self.__problem.write_solution(fn)

    #region Configuration
        
    def __get_target_class_config(self):
        target_classes = {}

        for name, options in self.__get_netflow_option('targetClasses', {}).items():
            
            min_targets = options['min_targets'] if 'min_targets' in options else None
            max_targets = options['max_targets'] if 'max_targets' in options else None
            non_observation_cost = options['non_observation_cost'] if 'non_observation_cost' in options else None
            partial_observation_cost = options['partial_observation_cost'] if 'partial_observation_cost' in options else None
            calib = options['calib'] if 'calib' in options else False
            
            # sanity check for science targets: make sure that partialObservationCost
            # is larger than nonObservationCost
            if not calib \
                and non_observation_cost is not None \
                and partial_observation_cost is not None \
                and partial_observation_cost < non_observation_cost:
                
                raise ValueError(
                    f"Found target class `{name}` where partialObservationCost "
                    "is smaller than nonObservationCost")

            target_classes[name] = SimpleNamespace(
                min_targets = min_targets,
                max_targets = max_targets,
                non_observation_cost = non_observation_cost,
                partial_observation_cost = partial_observation_cost,
                calib = calib,
            )

        return target_classes
    
    def __get_forbidden_pairs_config(self):
        forbidden_pairs = self.__get_netflow_option('forbiddenPairs', None)
        fpp = []
        if forbidden_pairs is not None:
            # Reindex targets for faster search
            df = self.__targets.set_index('id')

            for i, pair in enumerate(forbidden_pairs):
                if len(pair) == 0 or len(pair) > 2:
                    raise ValueError(f"Found an incorrect number of target ids in forbidden pair list at index {i}.")
            
                tidx_list = [ df.index.get_loc(p) for p in pair ]
                fpp.append(tidx_list)
                
        return fpp

    def __get_black_dot_config(self):
        black_dot_penalty = self.__get_netflow_option('blackDotPenalty', None)
        
        if black_dot_penalty is not None:
            black_dots = SimpleNamespace(
                # Closest black dots for each cobra
                black_dot_list = self.__positioner.nf_get_closest_dots(),
                black_dot_penalty = black_dot_penalty
            )
        else:
            black_dots = None

        return black_dots
    
    def __get_cobra_groups_config(self):
        cobra_groups = {}
        
        for name, options in self.__get_netflow_option('cobraGroups', {}).items():
            for item in ['groups', 'target_classes', 'non_observation_cost']:
                if item not in options:
                    raise RuntimeError(f'Config entry `{item}` is missing for cobra group `{name}`.')
                
            groups = options['groups']
            max_group = np.max(groups)

            min_targets = options['min_targets'] if 'min_targets' in options else None
            max_targets = options['max_targets'] if 'max_targets' in options else None

            if min_targets is None and max_targets is None:
                raise RuntimeError(f'Config entry `min_targets` and `max_targets` are missing for cobra group `{name}`.')

            cobra_groups[name] = SimpleNamespace(
                groups = groups,
                max_group = max_group,
                target_classes = set(options['target_classes']),
                min_targets = min_targets,
                max_targets = max_targets,
                non_observation_cost = options['non_observation_cost'],
            )

        return cobra_groups
    
    def __get_time_budget_config(self):
        time_budgets = {}

        for name, options in self.__get_netflow_option('timeBudgets', {}).items():
            for item in ['target_classes', 'budget']:
                if item not in options:
                    raise RuntimeError(f'Config entry `{item}` is missing for time budget `{name}`.')
                
            time_budgets[name] = SimpleNamespace(
                target_classes = set(options['target_classes']),
                budget = options['budget'],
            )

        return time_budgets
    
    #endregion
        
    def __append_targets(self, catalog, id_column, prefix, exp_time=None, priority=None, penalty=None, mask=None, filter=None):
        df = catalog.get_data(mask=mask, filter=filter)

        # TODO: add pm, epoch
        targets = pd.DataFrame({ 'id': df[id_column],
                                 'RA': df['RA'],
                                 'Dec': df['Dec'] })
        
        # This is a per-object penalty for observing calibration targets
        if penalty is not None:
            targets['penalty'] = penalty
        elif 'penalty' in df.columns:
            targets['penalty'] = df['penalty']
        else:
            targets['penalty'] = 0

        targets['prefix'] = prefix

        if prefix == 'sci':
            if exp_time is not None:
                targets['exp_time'] = exp_time
            else:
                targets['exp_time'] = df['exp_time']

            if priority is not None:
                targets['priority'] = priority
            else:
                targets['priority'] = df['priority']

            targets['class'] = targets[['prefix', 'priority']].apply(lambda r: f"{r['prefix']}_P{r['priority']}", axis=1)
        else:
            targets['class'] = prefix

            # Calibration targets have no prescribed exposure time and priority
            targets['exp_time'] = -1
            targets['priority'] = -1

        # Append to the existing list
        if self.__targets is None:
            self.__targets = targets
        else:
            self.__targets = pd.concat([self.__targets, targets])
            
    def append_science_targets(self, catalog, exp_time=None, priority=None, mask=None, filter=None,):
        """Add science targets"""

        self.__append_targets(catalog, 'objid', 'sci', exp_time=exp_time, priority=priority, mask=mask, filter=filter)

    def append_sky_targets(self, sky, mask=None, filter=None):
        """Add sky positions"""

        self.__append_targets(sky, 'skyid', prefix='sky', mask=mask, filter=filter)

    def append_fluxstd_targets(self, fluxstd, mask=None, filter=None):
        """Add flux standard positions"""

        self.__append_targets(fluxstd, 'objid', prefix='cal', mask=mask, filter=filter)

    def build(self):
        """Construct the ILP problem"""
        # Load configuration
        logging.info("Processing configuration")

        self.__positioner = self.__positioner_type()

        # Target classes
        self.__target_classes = self.__get_target_class_config()

        # Forbidden target pairs
        self.__forbidden_pairs = self.__get_forbidden_pairs_config()

        # Cobras positioned too close to a black dot can get a penalty
        self.__black_dots = self.__get_black_dot_config()

        # Cobra groups are defined to set a minimum number of calibration targets in each
        self.__cobra_groups = self.__get_cobra_groups_config()

        # Science program time budgets
        self.__time_budgets = self.__get_time_budget_config()

        # Optimize data access
        self.__calculate_exp_time()
        self.__calculate_target_visits()
        self.__cache_targets()

        # Build the problem
        self.__calculate_target_fp_pos()
        self.__build_ilp_problem()
                    
    def __calculate_exp_time(self):
        # Since targeting is done in fix quanta in time, make sure all pointings and visits
        # have the same exposure time

        self.__visit_exp_time = None
        for pointing in self.__visits:
            if self.__visit_exp_time is None:
                self.__visit_exp_time = pointing.exp_time
            elif self.__visit_exp_time != pointing.exp_time:
                raise RuntimeError("Exposure time of every pointing and visit should be the same.")
            
    def __calculate_target_visits(self):
        """
        Calculate the number of required visits for each target.

        Assumes the same exposure time for each visit.
        """

        sci = (self.__targets['prefix'] == 'sci')

        # Reset the number of done visits for non science targets
        if 'done_visits' in self.__targets.columns:
            self.__targets['done_visits'][~sci]
        else:
            self.__targets['done_visits'] = 0

        # Calculate the number of required visits from the exposure time
        # The default is 0 for the non-science targets
        self.__targets['req_visits'] = 0
        self.__targets['req_visits'][sci] = np.int64(np.ceil(self.__targets['exp_time'][sci] / self.__visit_exp_time))

        pass

    def __cache_targets(self):
        # Extract the contents of the Pandas DataFrame for faster indexed access
        
        # TODO: add pm, epoch, etc. required by coordinate transform

        self.__target_cache = SimpleNamespace(
            id = np.array(self.__targets['id'].astype(np.int64)),
            ra = np.array(self.__targets['RA'].astype(np.float64)),
            dec = np.array(self.__targets['Dec'].astype(np.float64)),
            target_class = np.array(self.__targets['class'].astype(str)),
            prefix = np.array(self.__targets['prefix'].astype(str)),
            req_visits = np.array(self.__targets['req_visits'].astype(np.int32)),
            done_visits = np.array(self.__targets['done_visits'].astype(np.int32)),
            penalty = np.array(self.__targets['penalty'].astype(np.int32)),
        )

    def __calculate_target_fp_pos(self):
        """
        Calculate focal plane positions, etc. for each visit
        """

        telescopes = []
        self.__target_fp_pos = []
        for i, pointing in enumerate(self.__visits):
            # Telescope of the netflow lib is equivalent of Pointing + PFI of this lib
            telescope = self.__telescope_type(pointing)            
            telescopes.append(telescope)

            # Calculate focal plane positions of targets
            fpp, _ = telescope.world_to_fp_pos(self.__target_cache.ra, self.__target_cache.dec)
            self.__target_fp_pos.append(fpp[..., 0] + 1j * fpp[..., 1])

    def __make_name(self, *parts):
        ### SLOW ### 4M calls!
        return "_".join(str(x) for x in parts)

        # name = hex(self.__name_counter)
        # self.__name_counter += 1
        # return name

    def __add_variable(self, name, lo=None, hi=None):
        v = self.__problem.add_variable(name, lo=lo, hi=hi)
        self.__variables.all[name] = v
        return v

    def __add_cost(self, cost):
        self.__cost.append(cost)
        self.__problem.add_cost(cost)

    def __add_constraint(self, name, constraint):
        self.__constraints.all[name] = constraint
        self.__problem.add_constraint(name, constraint)

    def __build_ilp_problem(self):
        """
        Construct the ILP problem by defining the variables and constraints.
        """

        # Number of visits
        nvisits = len(self.__visits)

        self.__problem = self.__problem_type(self.__name, self.__solver_options)

        self.__variables = SimpleNamespace(
            all = dict(),

            CG_o = defaultdict(list),        # Cobra group outflows, key: (cobra_group_name, visit_idx, cobra_group)
            CG_Tv_Cv = defaultdict(list),    # Tv_Cv variables relevant to each cobra group, for each visit

            Cv_i = defaultdict(list),        # Cobra visit inflows, key: (cobra_idx, visit_idx)
            CTCv_o = defaultdict(list),      # Calibration target class visit outflows, key: (target_class, visit_idx)
            T_o = defaultdict(list),         # Target outflows (only science targets)
            T_i = defaultdict(list),         # Target inflows (only science targets), key: (target_idx)
            Tv_o = defaultdict(list),        # Target visit outflows, key: (target_idx, visit_idx)
            Tv_i = defaultdict(list),        # Target visit inflows, key: (target_idx, visit_idx)
            STC_o = defaultdict(list),       # Science Target outflows, key: (target_class)
        )

        self.__constraints = SimpleNamespace(
            all = dict(),

            Tv_o_coll = dict(),             # Fiber (endpoint or elbow) collision constraints,
                                            #      key: (tidx1, tidx2, visit_idx) if endpoint collisions
                                            #      key: (cidx1, cidx2, visit_idx) if elbow collisions -- why?
            Tv_o_forb = dict(),             # Forbidden pairs, key: (tidx1, tidx2, visit_idx)

            CG_Tv_Cv_min = dict(),          # Cobra group minimum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            CG_Tv_Cv_max = dict(),          # Cobra group maximum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            Cv_i_sum = dict(),              # At most one target per cobra, key: (cobra_idx, visit_idx)
            Cv_i_max = dict(),              # Hard upper limit on the number of assigned fibers, key (visit_idx)
            Tv_i_Tv_o_sum = dict(),         # Target visit nodes must be balanced, key: (target_id, visit_idx)
            T_i_T_o_sum = dict(),           # Inflow and outflow at every T node must be balanced, key: (target_id)
            CTCv_o_min = dict(),            # At least a required number of calibration target in each class (target), key: (target_class, visit_idx)
            STC_o_sum = dict(),             # Max number of science targets per target class, kez: (target_class)
            Tv_i_sum = dict(),              # Science program time budget, key: (budget_idx)
        )

        self.__cost = []

        logging.info("Creating network topology")

        # Create science target variables that are independent of visit
        self.__create_science_target_class_variables()
        self.__create_science_target_variables()
        
        # Create the variables for each visit
        for ivis in range(nvisits):
            # Create cobra group sinks for each cobra group within each
            # cobra group definition
            self.__create_cobra_group_variables(ivis)

            # Create calibration target class sinks and define cost
            # for each calibration target class
            self.__create_calib_target_class_variables(ivis)

            logging.debug("Calculating visibilities")
            logging.debug(f"  exposure {ivis + 1}")
            vis_elbow = self.__positioner.nf_get_visibility_and_elbow(self.__target_fp_pos[ivis])
            self.__create_visit_variables(ivis, vis_elbow)
            self.__create_collision_constraints(ivis, vis_elbow)

            # Add constraints for forbidden pairs of targets
            if self.__forbidden_pairs is not None and len(self.__forbidden_pairs) > 0:
                logging.debug("Adding forbidden pair constraints")
                self.__create_forbidden_pairs_constraints(ivis)
        
        logging.info("Adding constraints")

        # TODO: group these up into functions?

        # Every cobra can observe at most one target per visit    
        self.__create_Cv_i_constraints()

        # Inflow and outflow at every Tv node must be balanced
        self.__create_Tv_i_Tv_o_constraints()

        # Inflow and outflow at every Tv node must be balanced
        self.__create_T_i_T_o_constraints()
            
        # Every calibration target class must be observed a minimum number of times
        # every visit
        self.__create_CTCv_o_constraints()

        # The maximum number of targets to be observed within a science target class
        # It defaults to the total number of targets but can be overridden for each target class
        # TODO: this could be configured in a more detailed manner
        self.__create_STC_o_constraints()

        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget
        self.__create_time_budget_constraints()

        # Make sure that there are enough sky targets in every Cobra location group
        # and instrument region
        self.__create_cobra_group_constraints(nvisits)

        # Make sure that enough fibers are kept unassigned, if this was requested
        self.__create_unassigned_fiber_constraints(nvisits)

    #region Variables and costs
        
    def __create_science_target_class_variables(self):
        # TODO: replace member `calib` with `prefix` to be consistent with targets
        for tckey, target_class in self.__target_classes.items():
            if not target_class.calib:
                # Science Target class node to sink
                self.__create_STC_sink(tckey)

    def __create_STC_sink(self, tckey):
        f = self.__add_variable(self.__make_name("STC_sink", tckey), 0, None)
        self.__variables.STC_o[tckey].append(f)
        self.__add_cost(f * self.__target_classes[tckey].non_observation_cost)
        
    def __create_science_target_variables(self):
        for tidx in range(len(self.__targets)):
            target_id = self.__target_cache.id[tidx]
            target_class = self.__target_cache.target_class[tidx]
            target_prefix = self.__target_cache.prefix[tidx]

            if target_prefix == 'sci':
                # Science Target class node to target node
                self.__create_STC_T(tidx, target_class, target_id)

                # Science Target node to sink
                self.__create_ST_sink(tidx, target_class, target_id)

    def __create_STC_T(self, tidx, target_class, target_id):
        # TODO: This is a bit fishy here. If the target is marked as already observed
        #       by done_visits > 0, why is a new variable added at all?
        lo = 1 if self.__target_cache.done_visits[tidx] > 0 else 0

        f = self.__add_variable(self.__make_name("STC_T", target_class, target_id), lo, 1)
        self.__variables.T_i[tidx].append(f)
        self.__variables.STC_o[target_class].append(f)

    def __create_ST_sink(self, tidx, target_class, target_id):
        # TODO: we could calculate a maximum for this, which is the total number of visits

        # TODO: this is wrong, a target is only partially observed if it doesn't have as
        #       many visits as required by observation time

        # TODO: need to have a constraint somewhere which makes sure "ST_sink"
        #       doesn't take values larger than 1
        #       partial observation cost applies only when the number of visits is less
        #       than required by the exposure time

        raise NotImplementedError()

        f = self.__add_variable(self.__make_name("ST_sink", target_id), 0, None)
        self.__variables.T_o[tidx].append(f)
        self.__add_cost(f * self.__target_classes[target_class].partial_observation_cost)   
            
    def __create_cobra_group_variables(self, ivis):
        # Cobra groups are defined by array with the same size as the number of cobras
        # Each item of this array is a cobra group number
        # There can be a minimum required number and maximum allowed number of targets
        # in each cobra group. This is to set a lower limit of sky fibers in every sky location
        # and instrument region.

        for name, options in self.__cobra_groups.items():
            # Create a sink for each group and each visit
            for cgidx in range(options.max_group + 1):
                # Add overflow arcs to the sink
                self.__create_CG_sink(ivis, cgidx, name, options)

    def __create_CG_sink(self, ivis, cgidx, name, options):
        # Add overflow arcs to the sink
        f = self.__add_variable(self.__make_name('CG_sink', name, cgidx, ivis), 0, None)
        self.__variables.CG_o[(name, cgidx, ivis)].append(f)
        self.__add_cost(f * options.non_observation_cost)

    def __create_calib_target_class_variables(self, ivis):
        # Calibration target class visit outflows, key: (target_class, visit_idx)
        # TODO: replace member `calib` with `prefix` to be consistent with targets
        for target_class, options in self.__target_classes.items():
            if options.calib:
                self.__create_CTCv_sink(ivis, target_class, options)
                
    def __create_CTCv_sink(self, ivis, target_class, options):
        f = self.__add_variable(self.__make_name("CTCv_sink", target_class, ivis), 0, None)
        self.__variables.CTCv_o[(target_class, ivis)].append(f)
        self.__add_cost(f * options.non_observation_cost)
    
    def __create_visit_variables(self, ivis, vis_elbow):
        # For each target-cobra pair with elbow position
        for tidx, cidx_elbow in vis_elbow.items():
            target_id = self.__target_cache.id[tidx]
            target_class = self.__target_cache.target_class[tidx]
            target_prefix = self.__target_cache.prefix[tidx]

            if target_prefix == 'sci':
                # Target node to target visit node
                self.__create_T_Tv(ivis, tidx, target_id)
            elif target_prefix in ['sky', 'cal']:
                # Calibration Target class node to target visit node
                self.__create_CTCv_Tv(ivis, tidx, target_class, target_id)
            else:
                raise NotImplementedError()

            for (cidx, _) in cidx_elbow:
                # Target visit node to cobra visit node
                self.__create_Tv_Cv(ivis, tidx, cidx, target_class, target_prefix)

        # If requested, penalize non-allocated fibers
        # Sum up the all Cv_i edges for the current visit and penalize its difference from
        # the total number of cobras
        fiber_non_allocation_cost = self.__get_netflow_option('fiberNonAllocationCost', 0)
        if fiber_non_allocation_cost != 0:
            # TODO: consider storing Cv_i organized by visit instead of cobra as well
            relevant_vars = [ var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == ivis ]
            relevant_vars = [ item for sublist in relevant_vars for item in sublist ]
            self.__add_cost(fiber_non_allocation_cost *
                            (self.__positioner.bench.cobras.nCobras - self.__problem.sum(relevant_vars)))

    def __create_T_Tv(self, ivis, tidx, target_id):
        # Target node to target visit node
        f = self.__add_variable(self.__make_name("T_Tv", target_id, ivis), 0, 1)
        self.__variables.T_o[tidx].append(f)
        self.__variables.Tv_i[(tidx, ivis)].append(f)

    def __create_CTCv_Tv(self, ivis, tidx, target_class, target_id):
        # Calibration Target class node to target visit node
        f = self.__add_variable(self.__make_name("CTCv_Tv", target_class, target_id, ivis), 0, 1)
        self.__variables.Tv_i[(tidx, ivis)].append(f)
        self.__variables.CTCv_o[(target_class, ivis)].append(f)

        # TODO: Why do we penalize observing a particular calibration target?
        if self.__target_cache.penalty[tidx] != 0:
            self.__add_cost(f * self.__target_cache.penalty[tidx])

    def __create_Tv_Cv(self, ivis, tidx, cidx, target_class, target_prefix):
        f = self.__add_variable(self.__make_name("Tv_Cv", tidx, cidx, ivis), 0, 1)
        self.__variables.Cv_i[(cidx, ivis)].append(f)
        self.__variables.Tv_o[(tidx, ivis)].append((f, cidx))

        # Save the variable to the list of each cobra group to which it's relevant
        for cg_name, options in self.__cobra_groups.items():
            if target_class in options.target_classes:
                self.__variables.CG_Tv_Cv[(cg_name, ivis, options.groups[cidx])].append(f)

        # Cost of the visit
        total_cost = self.__visits[ivis].obs_cost

        # Cost of moving the cobra
        cobra_move_cost = self.__get_netflow_option('cobraMoveCost', None)
        if cobra_move_cost is not None:
            dist = np.abs(self.__positioner.bench.cobras.centers[cidx] - self.__target_fp_pos[ivis][tidx])
            total_cost += cobra_move_cost(dist)
        
        # Cost of closest black dots for each cobra
        if self.__black_dots is not None:
            dist = np.min(np.abs(self.__black_dots.black_dot_list[cidx] - self.__target_fp_pos[ivis][tidx]))
            total_cost += self.__black_dots.black_dot_penalty(dist)
        
        if total_cost != 0:
            self.__add_cost(f * total_cost)

    #endregion
    #region Constraints

    def __create_collision_constraints(self, ivis, vis_elbow):
        # TODO: create individual function for each type of constraint

        # Add constraints 
        logging.debug(f"Adding constraints for visit {ivis}")

        # Avoid endpoint or elbow collisions
        collision_distance = self.__get_netflow_option('collision_distance', 0.0)
        elbow_collisions = self.__get_netflow_option('elbow_collisions', False)
        
        if collision_distance > 0.0:
            if not elbow_collisions:
                logging.debug("Adding endpoint collision constraints")
                self.__create_endpoint_collision_constraints(ivis, vis_elbow, collision_distance)
            else:
                logging.debug("Adding elbow collision constraints")
                self.__create_elbow_collision_constraints(ivis, vis_elbow, collision_distance)
                
    def __create_endpoint_collision_constraints(self, ivis, vis_elbow, collision_distance):
        # Determine target indices visible by this cobra and its neighbors
        ignore_endpoint_collisions = self.__get_debug_option('ignoreEndpointCollisions', False)

        colliding_pairs = self.__positioner.nf_get_colliding_pairs(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
        keys = self.__variables.Tv_o.keys()
        keys = set(key[0] for key in keys if key[1] == ivis)
        for p in colliding_pairs:
            if p[0] in keys and p[1] in keys:
                flows = [ v for v, cidx in self.__variables.Tv_o[(p[0], ivis)] ] + \
                        [ v for v, cidx in self.__variables.Tv_o[(p[1], ivis)] ]
                name = self.__make_name("Tv_o_coll", p[0], p[1], ivis)
                constr = self.__problem.sum(flows) <= 1
                self.__constraints.Tv_o_coll[(p[0], p[1], ivis)] = constr
                if not ignore_endpoint_collisions:
                    self.__add_constraint(name, constr)

    def __create_elbow_collision_constraints(self, ivis, vis_elbow, collision_distance):
        # Determine targets accessible by two different cobras that would cause an elbow collision
        # colliding_elbows contains a list of tidx keyed by cidx and tidx

        ignore_elbow_collisions = self.__get_debug_option('ignoreElbowCollisions', False)

        ### SLOW ###

        colliding_elbows = self.__positioner.nf_get_colliding_elbows(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
        for (cidx1, tidx1), tidx2_list in colliding_elbows.items():
            for f, cidx2 in self.__variables.Tv_o[(tidx1, ivis)]:
                if cidx2 == cidx1:
                    flow0 = f

            for tidx2 in tidx2_list:
                if True:  # idx2 != tidx1:
                    flow = [ flow0 ]
                    flow += [ f for f, cidx2 in self.__variables.Tv_o[(tidx2, ivis)] if cidx2 != cidx1 ]
        
                    name = self.__make_name("Tv_o_coll", tidx1, cidx1, tidx2, cidx2, ivis)
                    constr = self.__problem.sum(flow) <= 1
                    self.__constraints.Tv_o_coll[(tidx1, cidx1, tidx2, cidx2, ivis)] = constr
                    if not ignore_elbow_collisions:
                        self.__add_constraint(name, constr)

    def __create_forbidden_pairs_constraints(self, ivis):
        # All edges of visible targets, relevant for this visit

        ignore_forbidden_pairs = self.__get_debug_option('ignoreForbiddenPairs', False)
        ignore_forbidden_singles = self.__get_debug_option('ignoreForbiddenSingles', False)

        tidx_set = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == ivis)
        for p in self.__forbidden_pairs:
            if len(p) == 2:
                [tidx1, tidx2] = p
                if tidx1 in tidx_set and tidx2 in tidx_set:
                    flows = [ v for v, cidx in self.__variables.Tv_o[(tidx1, ivis)] ] + \
                            [ v for v, cidx in self.__variables.Tv_o[(tidx2, ivis)] ]
                    name = self.__make_name("Tv_o_forb", tidx1, tidx2, ivis)
                    constr = self.__problem.sum(flows) <= 1
                    self.__constraints.Tv_o_forb[(tidx1, tidx2, ivis)] = constr
                    if not ignore_forbidden_pairs:
                        self.__add_constraint(name, constr)
            elif len(p) == 1:
                [tidx] = p
                if tidx in tidx_set:
                    flows = [ v for v, cidx in self.__variables.Tv_o[(tidx, ivis)] ]
                    name = self.__make_name("Tv_o_forb", tidx, tidx, ivis)
                    constr = self.__problem.sum(flows) == 0
                    self.__constraints.Tv_o_forb[(tidx, tidx, ivis)] = constr
                    if not ignore_forbidden_singles:
                        self.__add_constraint(name, constr)
            else:
                raise RuntimeError("oops")

    def __create_Cv_i_constraints(self):
        # Every cobra can observe at most one target per visit
        for (cidx, vidx), inflow in self.__variables.Cv_i.items():
            name = self.__make_name("Cv_i_sum", cidx, vidx)
            constr = self.__problem.sum([ f for f in inflow ]) <= 1
            self.__constraints.Cv_i_sum[(cidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_Tv_i_Tv_o_constraints(self):
        # Inflow and outflow at every Tv node must be balanced
        for (tidx, vidx), ivars in self.__variables.Tv_i.items():
            ovars = self.__variables.Tv_o[(tidx, vidx)]
            target_id = self.__target_cache.id[tidx]
            name = self.__make_name("Tv_i_Tv_o_sum", target_id, vidx)
            constr = self.__problem.sum([ v for v in ivars ] + [ -v[0] for v in ovars ]) == 0
            self.__constraints.Tv_i_Tv_o_sum[(tidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_T_i_T_o_constraints(self):
        # Inflow and outflow at every T node must be balanced
        for tidx, ivars in self.__variables.T_i.items():
            ovars = self.__variables.T_o[tidx]
            target_id = self.__target_cache.id[tidx]
            nvis = max(0, self.__target_cache.req_visits[tidx] - self.__target_cache.done_visits[tidx])
            name = self.__make_name("T_i_T_o_sum", target_id)
            constr = self.__problem.sum([ nvis * v for v in ivars ] + [ -v for v in ovars ]) == 0
            self.__constraints.T_i_T_o_sum[(target_id)] = constr
            self.__add_constraint(name, constr)

    def __create_CTCv_o_constraints(self):
        # Every calibration target class must be observed a minimum number of times
        # every visit

        ignore_calib_target_class_minimum = self.__get_debug_option('ignoreCalibTargetClassMinimum', False)

        for (target_class, vidx), vars in self.__variables.CTCv_o.items():
            min_targets = self.__target_classes[target_class].min_targets
            if min_targets is not None:
                name = self.__make_name("CTCv_o_min", target_class, vidx)
                constr = self.__problem.sum([ v for v in vars ]) >= min_targets
                self.__constraints.CTCv_o_min[(target_class, vidx)] = constr
                if not ignore_calib_target_class_minimum:
                    self.__add_constraint(name, constr)

    def __create_STC_o_constraints(self):
        # Science targets must be either observed or go to the sink

        ignore_science_target_class_total = self.__get_debug_option('ignoreScienceTargetClassTotal', False)

        for target_class, vars in self.__variables.STC_o.items():
            if self.__target_classes[target_class].max_targets is None:
                # STC_o contains variables for each target in the science target class plus
                # one that goes to the sink
                max_targets = len(vars) - 1
            else:
                # TODO: Is this correct here? I think the number of science targets
                #       cannot be limited this way
                raise NotImplementedError()
                max_targets = self.__target_classes[target_class].max_targets
                
            name = self.__make_name("STC_o_sum", target_class)
            constr = self.__problem.sum([ v for v in vars ]) == max_targets
            self.__constraints.STC_o_sum[target_class] = constr
            if not ignore_science_target_class_total:
                self.__add_constraint(name, constr)

    def __create_time_budget_constraints(self):
        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget

        ignore_time_budget = self.__get_debug_option('ignoreTimeBudget', False)

        for budget_name, options in self.__time_budgets.items():
            budget_variables = []
            # Collect all visits of targets that fall into to budget
            for (tidx, ivis), val in self.__variables.Tv_i.items():
                target_class = self.__target_cache.target_class[tidx]
                if target_class in options.target_classes:
                    budget_variables.append(val[0])     # TODO: Why the index? Apparently a single var is stored in a list

            # TODO: This assumes a single per visit exposure time which is fine for now
            #       but the logic could be extended further
            name = self.__make_name("Tv_i_sum", budget_name)
            constr = self.__visit_exp_time * self.__problem.sum([ v for v in budget_variables ]) <= 3600 * options.budget
            self.__constraints.Tv_i_sum[budget_name] = constr
            if not ignore_time_budget:
                self.__add_constraint(name, constr)

    def __create_cobra_group_constraints(self, nvisits):
        # Make sure that there are enough targets in every cobra group for each visit

        ignore_cobra_group_minimum = self.__get_debug_option('ignoreCobraGroupMinimum', False)
        ignore_cobra_group_maximum = self.__get_debug_option('ignoreCobraGroupMaximum', False)

        for name, options in self.__cobra_groups.items():
            for ivis in range(nvisits):
                for i in range(options.max_group + 1):
                    variables = self.__variables.CG_Tv_Cv[(name, ivis, i)]
                    if len(variables) > 0:
                        if options.min_targets is not None:
                            name = self.__make_name("CG_Tv_Cv_min", name, ivis, i)
                            constr = self.__problem.sum([ v for v in variables ]) >= options.min_targets
                            self.__constraints.CG_Tv_Cv_min[name, ivis, i] = constr
                            if not ignore_cobra_group_minimum:
                                self.__add_constraint(name, constr)
                        if options.max_targets is not None:
                            name = self.__make_name("CG_Tv_Cv_max", name, ivis, i)
                            constr = self.__problem.sum([ v for v in variables ]) <= options.min_targets
                            self.__constraints.CG_Tv_Cv_max[name, ivis, i] = constr
                            if not ignore_cobra_group_maximum:
                                self.__add_constraint(name, constr)

    def __create_unassigned_fiber_constraints(self, nvisits):
        # Make sure that enough fibers are kept unassigned, if this was requested
        # This is done by setting an upper limit on the sum of Cv_i edges

        ignore_unassigned_minimum = self.__get_debug_option('ignoreUnassignedMinimum', False)
        numReservedFibers = self.__get_netflow_option('numReservedFibers', 0)

        if numReservedFibers > 0:
            maxAssignableFibers = self.__positioner.bench.cobras.nCobras - numReservedFibers
            for ivis in range(nvisits):
                variables = [var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == ivis]
                variables = [item for sublist in variables for item in sublist]
                name = self.__make_name("Cv_i_max", ivis)
                constr = self.__problem.sum(variables) <= maxAssignableFibers
                self.__constraints.Cv_i_max[ivis] = constr
                if not ignore_unassigned_minimum:
                    self.__add_constraint(name, constr)

    #endregion

    def solve(self):
        self.__problem.solve()
        self.__extract_assignments()

    def __extract_assignments(self):
        """Extract the fiber assignments from an LP solution"""

        self.__assignments = [ {} for _ in range(len(self.__visits)) ]

        for k1, v1 in self.__problem._variables.items():
            if k1.startswith("Tv_Cv_"):
                visited = self.__problem.get_value(v1) > 0
                if visited:
                    _, _, tidx, cidx, ivis = k1.split("_")
                    self.__assignments[int(ivis)][int(tidx)] = int(cidx)

    def __get_idcol(self, catalog):
        for c in ['objid', 'skyid', 'calid']:
            if c in catalog.data.columns:
                return c
        
        raise RuntimeError()
    
    def get_fiber_assignments(self, catalog):
        """Returns a mask that indexes the selected targets within a catalog."""

        idcol = self.__get_idcol(catalog)

        assignments = []
        for i, p in enumerate(self.__visits):
            idx = np.array([ k for k in self.__assignments[i] ])                   # Target indices

            targets = pd.DataFrame(
                    {
                        f'target_{idcol}': [ np.int64(self.__targets.iloc[ti]['id']) for ti in idx ],
                        'fiberid': np.array([ self.__assignments[i][k] for k in idx ]) 
                    }
                ).astype({ 'fiberid': pd.Int32Dtype() }).set_index(f'target_{idcol}')

            # Find the index of each object in the catalog and return the matching fiberid
            assign = catalog.data[[idcol]].set_index(idcol).join(targets, how='left').reset_index()
            assignments.append(assign)

        return assignments

    def get_fiber_assignments_masks(self, catalog):
        idcol = self.__get_idcol(catalog)
        assignments = self.get_fiber_assignments(catalog)
        masks = []
        fiberids = []
        for i, p in enumerate(self.__visits):
            # Boolean array of targets that have fibers assigned to
            # Join is required because assignments might be ordered differently than catalog.data
            a = catalog.data[[idcol]].join(assignments[i].set_index(idcol), on=idcol)
            m = ~pd.isna(a['fiberid'])
            masks.append(np.array(m))
            fiberids.append(np.array(np.int32(a['fiberid'][m])))

        return masks, fiberids