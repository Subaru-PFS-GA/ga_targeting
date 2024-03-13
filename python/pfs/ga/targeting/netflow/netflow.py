import os
import logging
from typing import Callable
from collections import defaultdict
import numpy as np
import pandas as pd
from types import SimpleNamespace  

from ics.cobraOps.Bench import Bench

from ..util.args import *
from ..instrument import SubaruWFC, SubaruPFI
from .gurobiproblem import GurobiProblem

class Netflow():
    """
    Implements the Network Flow algorithm to optimize fiber allocation.
    """
    
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
        self.__bench = Bench(layout='full')

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
        
        self.__target_assignments = None
        self.__cobra_assignments = None
        self.__missed_targets = None
        self.__partially_observed_targets = None

    #region Property accessors

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

    def __get_target_assignments(self):
        return self.__target_assignments
    
    target_assignments = property(__get_target_assignments)

    def __get_cobra_assignments(self):
        return self.__cobra_assignments
    
    cobra_assignments = property(__get_cobra_assignments)

    #endregion

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
        
    #region Save and load

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

    #endregion
    #region Configuration
        
    def __get_target_class_config(self):
        target_classes = {}

        for name, options in self.__get_netflow_option('targetClasses', {}).items():
            
            prefix = options['prefix'] if 'prefix' in options else False
            min_targets = options['min_targets'] if 'min_targets' in options else None
            max_targets = options['max_targets'] if 'max_targets' in options else None
            non_observation_cost = options['non_observation_cost'] if 'non_observation_cost' in options else None
            partial_observation_cost = options['partial_observation_cost'] if 'partial_observation_cost' in options else None
            
            # Sanity check for science targets: make sure that partialObservationCost
            # is larger than nonObservationCost
            if prefix == 'sci' \
                and non_observation_cost is not None \
                and partial_observation_cost is not None \
                and partial_observation_cost < non_observation_cost:
                
                raise ValueError(
                    f"Found target class `{name}` where partialObservationCost "
                    "is smaller than nonObservationCost")

            target_classes[name] = SimpleNamespace(
                prefix = prefix,
                min_targets = min_targets,
                max_targets = max_targets,
                non_observation_cost = non_observation_cost,
                partial_observation_cost = partial_observation_cost,
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
                black_dot_list = self.__get_closest_black_dots(),
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
            ngroups = np.max(groups) + 1

            min_targets = options['min_targets'] if 'min_targets' in options else None
            max_targets = options['max_targets'] if 'max_targets' in options else None

            if min_targets is None and max_targets is None:
                raise RuntimeError(f'Config entry `min_targets` and `max_targets` are missing for cobra group `{name}`.')

            cobra_groups[name] = SimpleNamespace(
                groups = groups,
                ngroups = ngroups,
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
    #region PFI functions

    def __get_closest_black_dots(self):
        """
        For each cobra, return the list of focal plane black dot positions.
        
        Returns
        =======
        list : List of arrays of complex focal plane positions.
        """
        res = []
        for cidx in range(len(self.__bench.cobras.centers)):
            nb = self.__bench.getCobraNeighbors(cidx)
            res.append(self.__bench.blackDots.centers[[cidx] + list(nb)])

        return res
    
    def __get_visibility_and_elbow(self, fp_pos):
        """
        Calculate the visibility and the corresponding elbow position for each
        target of a single visit from the focal plane positions in ˙tpos˙.

        Parameters
        ==========

        tpos : numpy.ndarray
            Complex focal plane positions of each target.

        Returns
        =======
        dict : Dictionary of lists of (cobra ID, elbow position) pairs, keyed by target indices.
        """

        from ics.cobraOps.TargetGroup import TargetGroup
        from ics.cobraOps.TargetSelector import TargetSelector

        class DummyTargetSelector(TargetSelector):
            def run(self):
                return

            def selectTargets(self):
                return

        tgroup = TargetGroup(fp_pos)
        tselect = DummyTargetSelector(self.__bench, tgroup)
        tselect.calculateAccessibleTargets()
        targets = tselect.accessibleTargetIndices   # shape: (cobras, targets), padded with -1
        elbows = tselect.accessibleTargetElbows     # shape: (cobras, targets), padded with 0+0j
        
        # Build two dictionaries with different look-up directions
        targets_cobras = defaultdict(list)
        cobras_targets = defaultdict(list)

        for cidx in range(targets.shape[0]):
            for i, tidx in enumerate(targets[cidx, :]):
                if tidx >= 0:
                    targets_cobras[tidx].append((cidx, elbows[cidx, i]))
                    cobras_targets[cidx].append((tidx, elbows[cidx, i]))

        return targets_cobras, cobras_targets
    
    def __get_colliding_pairs(self, fp_pos, vis_elbow, dist):
        """Return the list of target pairs that would cause fiber
        collision when observed by two neighboring fibers.
        
        Parameters
        ==========

        fp_pos : array
            Complex focal plane positions of the targets
        vis_elbow: dict
            Visibility, dict of (cobra ID, elbow position), keyed by target index.
        dist:
            Maximum distance causing a collision.

        Returns
        =======
        set : pairs of target indices that would cause fiber top collisions.
        """

        # Collect targets associated with each cobra, for each visit
        fp_pos = np.array(fp_pos)
        ivis = defaultdict(list)
        for tidx, thing in vis_elbow.items():
            for (cidx, _) in thing:
                ivis[cidx].append(tidx)

        pairs = set()
        for cidx, i1 in ivis.items():
            # Determine target indices visible by this cobra and its neighbors
            nb = self.__bench.getCobraNeighbors(cidx)
            i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
            i2 = np.concatenate((i1, i2))
            i2 = np.unique(i2).astype(int)
            d = np.abs(np.subtract.outer(fp_pos[i1], fp_pos[i2]))
            for m in range(d.shape[0]):
                for n in range(d.shape[1]):
                    if d[m][n] < dist:
                        if i1[m] < i2[n]:               # Only store pairs once
                            pairs.add((i1[m], i2[n]))
        return pairs
    
    def __get_colliding_elbows(self, fp_pos, vis_elbow, dist):
        """
        For each target-cobra pair, and the corresponding elbow position,
        return the list of other targets that are too close to the "upper arm" of
        the cobra.
        
        Parameters
        ==========

        fp_pos : array
            Complex focal plane positions of the targets
        vis_elbow: dict
            Visibility, dict of (cobra ID, elbow position), keyed by target index.
        dist:
            Maximum distance causing a collision.

        Returns
        =======
        dict : Dictionary of list of targets too close to the cobra indexed by
               all possible target-cobra pairs.
        """

        # TODO: speed up by vectorizing loops

        # vis contains the visible cobra indices and corresponding elbow positions
        # by target index. Invert this and build dictionaries indexed by cobra
        # indices that contain the lists of targets with corresponding elbow positions.
        ivis = defaultdict(list)
        epos = defaultdict(list)    # target_index, elbow position pairs keyed by cobra_index
        for tidx, cidx_elbow in vis_elbow.items():
            for (cidx, elbowpos) in cidx_elbow:
                ivis[cidx].append(tidx)
                epos[cidx].append((tidx, elbowpos))

        res = defaultdict(list)
        for cidx, tidx_elbow in epos.items():
            # Determine target indices visible by neighbors of this cobra
            nb = self.__bench.getCobraNeighbors(cidx)       # list of cobra_index
            tmp = [ epos[j] for j in nb if j in epos ]      # all targets visible by neighboring cobras
            if len(tmp) > 0:
                i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
                i2 = np.unique(i2).astype(int)

                # For each target visible by this cobra and the corresponding elbow
                # position, find all targets which are too close to the "upper arm"
                # of the cobra
                for tidx, elbowpos in tidx_elbow:
                    ebp = np.full(len(i2), elbowpos)
                    tp = np.full(len(i2), fp_pos[tidx])
                    ti2 = fp_pos[i2]
                    d = self.__bench.distancesToLineSegments(ti2, tp, ebp)
                    res[(cidx, tidx)] += list(i2[d < dist])

        return res

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

    def __cache_targets(self):
        """Extract the contents of the Pandas DataFrame for faster indexed access."""
        
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
    
    def __init_problem(self):
        self.__problem = self.__problem_type(self.__name, self.__solver_options)

    def __add_cost(self, cost):
        self.__problem.add_cost(cost)

    def __init_variables(self):
        self.__variables = SimpleNamespace(
            all = dict(),

            STC_o = defaultdict(list),       # Science target class outflows, key: (target_class)
            STC_sink = dict(),               # Science target class sink, key: (target_class)
            CTCv_o = defaultdict(list),      # Calibration target class visit outflows, key: (target_class, visit_idx)
            CTCv_sink = dict(),              # Calibration target class sink, key: (target_class, visit_idx)
            
            T_i = dict(),                    # Target inflows (only science targets), key: (target_idx)
            T_o = defaultdict(list),         # Target outflows (only science targets), key: (target_idx)
            T_sink = dict(),                 # Target sinks (only science targets) key: (target_idx)

            Tv_i = defaultdict(list),        # Target visit inflows, key: (target_idx, visit_idx)
            Tv_o = defaultdict(list),        # Target visit outflows, key: (target_idx, visit_idx)
            Cv_i = defaultdict(list),        # Cobra visit inflows, key: (cobra_idx, visit_idx)            

            CG_i = defaultdict(list),        # Cobra group inflow variables, for each visit, key: (name, ivis, gidx)
        )

    def __add_variable(self, name, lo=None, hi=None):
        v = self.__problem.add_variable(name, lo=lo, hi=hi)
        self.__variables.all[name] = v
        return v
    
    def __add_variable_array(self, name, indexes, lo=None, hi=None, cost=None):
        v = self.__problem.add_variable_array(name, indexes, lo=lo, hi=hi, cost=cost)
        self.__variables.all[name] = v
        return v
    
    def __init_constraints(self):
        self.__constraints = SimpleNamespace(
            all = dict(),

            STC_o_sum = dict(),             # Total number of science targets per target class, key: (target_class)
            STC_o_min = dict(),             # Min number of science targets per target class, key: (target_class)
            STC_o_max = dict(),             # Max number of science targets per target class, key: (target_class)

            CTCv_o_sum = dict(),            # Total number of calibration targets per target class, key: (target_class, vidx)
            CTCv_o_min = dict(),            # At least a required number of calibration targets in each class (target), key: (target_class, visit_idx)
            CTCv_o_max = dict(),            # Maximum number of calibration targets in each class, key: (target_class, visit_idx)

            Tv_o_coll = dict(),             # Fiber (endpoint or elbow) collision constraints,
                                            #      key: (tidx1, tidx2, visit_idx) if endpoint collisions
                                            #      key: (cidx1, cidx2, visit_idx) if elbow collisions -- why?
            Tv_o_forb = dict(),             # Forbidden pairs, key: (tidx1, tidx2, visit_idx)

            CG_min = dict(),                # Cobra group minimum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            CG_max = dict(),                # Cobra group maximum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            Cv_i_sum = dict(),              # At most one target per cobra, key: (cobra_idx, visit_idx)
            Cv_i_max = dict(),              # Hard upper limit on the number of assigned fibers, key (visit_idx)
            Tv_i_Tv_o_sum = dict(),         # Target visit nodes must be balanced, key: (target_id, visit_idx)
            T_i_T_o_sum = dict(),           # Inflow and outflow at every T node must be balanced, key: (target_id)
            
            Tv_i_sum = dict(),              # Science program time budget, key: (budget_idx)
        )

    def __add_constraint(self, name, constraint):
        self.__constraints.all[name] = constraint
        self.__problem.add_constraint(name, constraint)

    def __build_ilp_problem(self):
        """
        Construct the ILP problem by defining the variables and constraints.
        """

        # Number of visits
        nvisits = len(self.__visits)

        self.__init_problem()
        self.__init_variables()
        self.__init_constraints()

        logging.info("Creating network topology")

        # Create science target variables that are independent of visit because
        # we want multiple visits of the same target
        self.__create_science_target_class_variables()

        # Create science target variables, for each visit
        # > STC_T
        self.__create_science_target_variables()
        
        # Create the variables for each visit
        for ivis in range(nvisits):
            logging.info(f'Processing exposure {ivis + 1}.')

            # Create calibration target class variables and define cost
            # for each calibration target class. These are different for each visit
            # because we don't necessarily need the same calibration targets at each visit.
            # > CTCv_sink
            logging.debug("Creating calibration target class variables.")
            self.__create_calib_target_class_variables(ivis)

            logging.info("Calculating visibilities.")
            vis_targets_elbows, vis_cobras_targets = self.__get_visibility_and_elbow(self.__target_fp_pos[ivis])

            # > T_Tv, CTCv_Tv, Tv_Cv
            logging.debug("Creating target and cobra visit variables.")
            self.__create_visit_variables(ivis, vis_targets_elbows, vis_cobras_targets)

            logging.debug("Creating cobra collision constraints.")
            self.__create_cobra_collision_constraints(ivis, vis_targets_elbows)

            logging.debug("Adding cobra non-allocation cost terms.")
            self.__add_cobra_non_allocation_cost(ivis)

            # Add constraints for forbidden pairs of targets
            if self.__forbidden_pairs is not None and len(self.__forbidden_pairs) > 0:
                logging.debug("Adding forbidden pair constraints")
                self.__create_forbidden_pairs_constraints(ivis)
        
        logging.info("Adding constraints")

        # The maximum number of targets to be observed within a science target class
        self.__create_science_target_class_constraints()

        # Every calibration target class must be observed a minimum number of times
        # every visit
        self.__create_calibration_target_class_constraints()

        # Inflow and outflow at every Tv node must be balanced
        ##############
        self.__create_science_target_constraints()

        # Inflow and outflow at every Tv node must be balanced
        self.__create_Tv_i_Tv_o_constraints()

        # Every cobra can observe at most one target per visit    
        self.__create_Cv_i_constraints()

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
        for target_class in self.__target_classes.keys():
            if self.__target_classes[target_class].prefix == 'sci':
                # Science Target class node to sink
                self.__create_STC_sink(target_class)

    def __create_STC_sink(self, target_class):
        # Non-observed science target add to the cost with a coefficient defined
        # for the particular science target class. Higher priority targets should add
        # more to the cost to prefer targeting them.
        f = self.__add_variable(self.__make_name("STC_sink", target_class), 0, None)
        self.__variables.STC_sink[target_class] = f
        self.__add_cost(f * self.__target_classes[target_class].non_observation_cost)
        
    def __create_science_target_variables(self):
        # TODO: replace this with the logic to implement step-by-step targeting
        force_already_observed = self.__get_netflow_option('forceAlreadyObserved', False)

        for tidx in range(len(self.__targets)):
            target_id = self.__target_cache.id[tidx]
            target_class = self.__target_cache.target_class[tidx]
            target_prefix = self.__target_cache.prefix[tidx]
            target_penalty = self.__target_cache.penalty[tidx]
            
            if target_prefix == 'sci':
                # Science Target class node to target node
                self.__create_STC_T(tidx, target_class, target_id, target_penalty, force_already_observed)

                # Science Target node to sink
                self.__create_T_sink(tidx, target_class, target_id)

    def __create_STC_T(self, tidx, target_class, target_id, target_penalty, force_already_observed):
        # TODO: This is a bit fishy here. If the target is marked as already observed
        #       by done_visits > 0, why is a new variable added at all?
        #       It seems that this setting is trying to force observation of targets
        #       that have already been observed but this might cause problems, for
        #       example when there is a fiber collision.

        # If the STC_T edge is forced to be 1, it will eventually force the observation of
        # the particular target because the edges of T are always balanced.
        # The problem is that we don't know at which visits and by what fiber which is also
        # necessary to resume targeting.

        # TODO: just remove this and replace with logic to implement step-by-step targeting
        # if force_already_observed:
        #     raise NotImplementedError()
        #     lo = 1 if self.__target_cache.done_visits[tidx] > 0 else 0
        # else:
        #     lo = 0

        f = self.__add_variable(self.__make_name("STC_T", target_class, target_id), 0, 1)
        self.__variables.T_i[tidx] = f
        self.__variables.STC_o[target_class].append(f)

        # Cost of observing the science target
        if target_penalty != 0:
            self.__add_cost(f * target_penalty)

    def __create_T_sink(self, tidx, target_class, target_id):
        # TODO: we could calculate a maximum for this, which is the total number of visits

        # TODO: this is wrong, a target is only partially observed if it doesn't have as
        #       many visits as required by observation time

        # TODO: way to force the required number of observations only:
        #       - set constraints on the T_i or the sink
        #       - set a partial obs cost and a too many obs cost
        #         -- but how to set the cost to zero when within the bounds?
        # https://stackoverflow.com/questions/69904853/gurobipy-optimization-constraint-to-make-variable-value-to-be-greater-than-100
        # https://support.gurobi.com/hc/en-us/articles/4414392016529-How-do-I-model-conditional-statements-in-Gurobi

        f = self.__add_variable(self.__make_name("T_sink", target_id), 0, None)
        self.__variables.T_sink[tidx] = f
        self.__add_cost(f * self.__target_classes[target_class].partial_observation_cost)

    def __create_calib_target_class_variables(self, ivis):
        # Calibration target class visit outflows, key: (target_class, visit_idx)
        for target_class, options in self.__target_classes.items():
            if options.prefix in ['sky', 'cal']:
                self.__create_CTCv_sink(ivis, target_class, options)
                
    def __create_CTCv_sink(self, ivis, target_class, options):
        f = self.__add_variable(self.__make_name("CTCv_sink", target_class, ivis), 0, None)
        self.__variables.CTCv_sink[(target_class, ivis)] = f
        self.__add_cost(f * options.non_observation_cost)
    
    def __create_visit_variables(self, ivis, vis_targets_elbows, vis_cobras_targets):
        tidx = np.array(list(vis_targets_elbows.keys()), dtype=int)

        # Create science targets T_Tv edges in batch
        sci_mask = self.__target_cache.prefix[tidx] == 'sci'
        self.__create_T_Tv(ivis, tidx[sci_mask])

        # Create calibration targets CTCv_Tv edges
        cal_mask = (self.__target_cache.prefix[tidx] == 'sky') | (self.__target_cache.prefix[tidx] == 'cal')
        self.__create_CTCv_Tv(ivis, tidx[cal_mask])

        # For each Cv, generate the incoming Tv_Cv edges
        # These differ in number but can be vectorized for each cobra
        for cidx, tidx_elbow in vis_cobras_targets.items():
            tidx = np.array([ ti for ti, _ in tidx_elbow ], dtype=int)
            self.__create_Tv_Cv_CG(ivis, cidx, tidx)

    def __create_T_Tv(self, ivis, tidx):
        vars = self.__add_variable_array(self.__make_name('T_Tv', ivis), tidx, 0, 1)
        for ti in tidx:
            f = vars[ti]
            self.__variables.T_o[ti].append(f)
            self.__variables.Tv_i[(ti, ivis)].append(f)
            
    def __create_CTCv_Tv(self, ivis, tidx):
        cost = self.__target_cache.penalty[tidx]
        name = self.__make_name('CTCv_Tv', ivis)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]
            self.__variables.Tv_i[(ti, ivis)].append(f)
            self.__variables.CTCv_o[(target_class, ivis)].append(f)

    def __create_Tv_Cv_CG(self, ivis, cidx, tidx):
        # Calculate the cost for each target - cobra assignment
        cost = np.zeros_like(tidx, dtype=float)
        
        # Cost of a single visit
        if self.__visits[ivis].obs_cost is not None:
            cost += self.__visits[ivis].obs_cost

        # Cost of moving the cobra
        cobra_move_cost = self.__get_netflow_option('cobraMoveCost', None)
        if cobra_move_cost is not None:
            dist = np.abs(self.__bench.cobras.centers[cidx] - self.__target_fp_pos[ivis][tidx])
            cost += cobra_move_cost(dist)
        
        # Cost of closest black dots for each cobra
        if self.__black_dots is not None:
            dist = np.min(np.abs(self.__black_dots.black_dot_list[cidx] - self.__target_fp_pos[ivis][tidx]))
            cost += self.__black_dots.black_dot_penalty(dist)

        # Create LP variables
        name = self.__make_name("Tv_Cv", ivis, cidx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]

            self.__variables.Cv_i[(cidx, ivis)].append(f)
            self.__variables.Tv_o[(ti, ivis)].append((f, cidx))

            # Save the variable to the list of each cobra group to which it's relevant
            for cg_name, options in self.__cobra_groups.items():
                if target_class in options.target_classes:
                    self.__variables.CG_i[(cg_name, ivis, options.groups[cidx])].append(f)

    #endregion
    #region Special cost term
            
    def __add_cobra_non_allocation_cost(self, ivis):
        # If requested, penalize non-allocated fibers
        # Sum up the all Cv_i edges for the current visit and penalize its difference from
        # the total number of cobras
        cobra_non_allocation_cost = self.__get_netflow_option('fiberNonAllocationCost', 0)
        if cobra_non_allocation_cost != 0:
            # TODO: consider storing Cv_i organized by visit instead of cobra as well
            relevant_vars = [ var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == ivis ]
            relevant_vars = [ item for sublist in relevant_vars for item in sublist ]
            self.__add_cost(cobra_non_allocation_cost *
                            (self.__bench.cobras.nCobras - self.__problem.sum(relevant_vars)))
            
    #endregion
    #region Constraints

    def __create_cobra_collision_constraints(self, ivis, vis_elbow):
        # Add constraints 
        logging.debug(f"Adding constraints for visit {ivis}")

        # Avoid endpoint or elbow collisions
        collision_distance = self.__get_netflow_option('collisionDistance', 0.0)
        elbow_collisions = self.__get_netflow_option('elbowCollisions', False)
        
        if collision_distance > 0.0:
            if not elbow_collisions:
                logging.debug("Adding endpoint collision constraints")
                self.__create_endpoint_collision_constraints(ivis, vis_elbow, collision_distance)
            else:
                logging.debug("Adding elbow collision constraints")
                self.__create_elbow_collision_constraints(ivis, vis_elbow, collision_distance)
                
    def __create_endpoint_collision_constraints(self, ivis, vis_elbow, collision_distance):
        ignore_endpoint_collisions = self.__get_debug_option('ignoreEndpointCollisions', False)

        if not ignore_endpoint_collisions:
            # Determine target indices visible by this cobra and its neighbors
            colliding_pairs = self.__get_colliding_pairs(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
            keys = self.__variables.Tv_o.keys()
            keys = set(key[0] for key in keys if key[1] == ivis)
            for p in colliding_pairs:
                if p[0] in keys and p[1] in keys:
                    vars = [ v for v, cidx in self.__variables.Tv_o[(p[0], ivis)] ] + \
                           [ v for v, cidx in self.__variables.Tv_o[(p[1], ivis)] ]
                    name = self.__make_name("Tv_o_coll", p[0], p[1], ivis)
                    constr = self.__problem.sum(vars) <= 1
                    self.__constraints.Tv_o_coll[(p[0], p[1], ivis)] = constr
                    self.__add_constraint(name, constr)

    def __create_elbow_collision_constraints(self, ivis, vis_elbow, collision_distance):
        # Determine targets accessible by two different cobras that would cause an elbow collision
        # colliding_elbows contains a list of tidx keyed by cidx and tidx

        ignore_elbow_collisions = self.__get_debug_option('ignoreElbowCollisions', False)

        ### SLOW ###

        colliding_elbows = self.__get_colliding_elbows(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
        for (cidx1, tidx1), tidx2_list in colliding_elbows.items():
            for f, cidx2 in self.__variables.Tv_o[(tidx1, ivis)]:
                if cidx2 == cidx1:
                    var0 = f

            for tidx2 in tidx2_list:
                if True:  # idx2 != tidx1:
                    vars = [ var0 ]
                    vars += [ f for f, cidx2 in self.__variables.Tv_o[(tidx2, ivis)] if cidx2 != cidx1 ]
        
                    name = self.__make_name("Tv_o_coll", tidx1, cidx1, tidx2, cidx2, ivis)
                    constr = self.__problem.sum(vars) <= 1
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
            
    def __create_science_target_class_constraints(self):
        # Science targets must be either observed or go to the sink
        # If a maximum on the science target is set fo the target class, enforce that
        # in a separate constraint
        # It defaults to the total number of targets but can be overridden for each target class
        ignore_science_target_class_minimum = self.__get_debug_option('ignoreScienceTargetClassMinimum', False)
        ignore_science_target_class_maximum = self.__get_debug_option('ignoreScienceTargetClassMaximum', False)
        
        for target_class, vars in self.__variables.STC_o.items():
            sink = self.__variables.STC_sink[target_class]
            num_targets = len(vars)

            name = self.__make_name("STC_o_sum", target_class)
            constr = self.__problem.sum(vars + [ sink ]) == num_targets
            self.__constraints.STC_o_sum[target_class] = constr
            self.__add_constraint(name, constr)

            # If a minimum or maximum constraint is set on the number observed targets in this
            # class, enforce it through a maximum constraint on the sum of the outgoing edges not
            # including the sink
            min_targets = self.__target_classes[target_class].min_targets
            if not ignore_science_target_class_minimum and min_targets is not None:
                name = self.__make_name("STC_o_min", target_class)
                constr = self.__problem.sum(vars) >= min_targets
                self.__constraints.STC_o_min[target_class] = constr
                self.__add_constraint(name, constr)

            max_targets = self.__target_classes[target_class].max_targets
            if not ignore_science_target_class_maximum and max_targets is not None:
                name = self.__make_name("STC_o_max", target_class)
                constr = self.__problem.sum(vars) <= max_targets
                self.__constraints.STC_o_max[target_class] = constr
                self.__add_constraint(name, constr)

    def __create_calibration_target_class_constraints(self):
        # The sum of all outgoing edges must be equal the number of calibration targets within each class
        
        # Every calibration target class must be observed a minimum number of times every visit
        ignore_calib_target_class_minimum = self.__get_debug_option('ignoreCalibTargetClassMinimum', False)
        ignore_calib_target_class_maximum = self.__get_debug_option('ignoreCalibTargetClassMaximum', False)
        
        for (target_class, vidx), vars in self.__variables.CTCv_o.items():
            sink = self.__variables.CTCv_sink[(target_class, vidx)]
            num_targets = len(vars)

            name = self.__make_name("CTCv_o_sum", target_class, vidx)
            constr = self.__problem.sum(vars + [ sink ]) == num_targets
            self.__constraints.CTCv_o_sum[(target_class, vidx)] = constr
            self.__add_constraint(name, constr)

            min_targets = self.__target_classes[target_class].min_targets
            if not ignore_calib_target_class_minimum and min_targets is not None:
                name = self.__make_name("CTCv_o_min", target_class, vidx)
                constr = self.__problem.sum(vars) >= min_targets
                self.__constraints.CTCv_o_min[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)

            max_targets = self.__target_classes[target_class].max_targets
            if not ignore_calib_target_class_maximum and max_targets is not None:
                name = self.__make_name("CTCv_o_max", target_class, vidx)
                constr = self.__problem.sum(vars) <= max_targets
                self.__constraints.CTCv_o_max[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)
            
    def __create_science_target_constraints(self):
        # Inflow and outflow at every T node must be balanced
        for tidx, in_flow in self.__variables.T_i.items():
            sink = self.__variables.T_sink[tidx]
            out_flow = self.__variables.T_o[tidx]
            target_id = self.__target_cache.id[tidx]
            req_visits = self.__target_cache.req_visits[tidx]

            name = self.__make_name("T_i_T_o_sum", target_id)
            constr = self.__problem.sum([ req_visits * in_flow ] + [ -v for v in out_flow ] + [ sink ]) == 0
            self.__constraints.T_i_T_o_sum[(target_id)] = constr
            self.__add_constraint(name, constr)

        # TODO: delete, once rewritten to step-by-step execution
        # constrain_already_observed = self.__get_netflow_option('constrainAlreadyObserved', False)
        # if constrain_already_observed:
        #     # TODO: replace this with logic to set a few Tv and Cv nodes to a fixed value to support
        #     #       iterative fitting

        #     # TODO: this enforces the number of required visits, only use this when forceAlreadyObserved,
        #     #       otherwise just add the cost terms
        #     raise NotImplementedError()

        #     # Inflow and outflow at every T node must be balanced
        #     for tidx, in_vars in self.__variables.T_i.items():
        #         out_vars = self.__variables.T_o[tidx]
        #         target_id = self.__target_cache.id[tidx]
        #         nvis = max(0, self.__target_cache.req_visits[tidx] - self.__target_cache.done_visits[tidx])
        #         name = self.__make_name("T_i_T_o_sum", target_id)
        #         constr = self.__problem.sum([ nvis * v for v in in_vars ] + [ -v for v in out_vars ]) == 0
        #         self.__constraints.T_i_T_o_sum[(target_id)] = constr
        #         self.__add_constraint(name, constr)
        # else:
        #     pass

    def __create_Tv_i_Tv_o_constraints(self):
        # Inflow and outflow at every Tv node must be balanced
        for (tidx, vidx), in_vars in self.__variables.Tv_i.items():
            out_vars = self.__variables.Tv_o[(tidx, vidx)]
            target_id = self.__target_cache.id[tidx]

            name = self.__make_name("Tv_i_Tv_o_sum", target_id, vidx)
            # TODO: why the index in outvars?
            constr = self.__problem.sum([ v for v in in_vars ] + [ -v[0] for v in out_vars ]) == 0
            self.__constraints.Tv_i_Tv_o_sum[(tidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_Cv_i_constraints(self):
        # Every cobra can observe at most one target per visit
        for (cidx, vidx), inflow in self.__variables.Cv_i.items():
            name = self.__make_name("Cv_i_sum", cidx, vidx)
            constr = self.__problem.sum([ f for f in inflow ]) <= 1
            self.__constraints.Cv_i_sum[(cidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_time_budget_constraints(self):
        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget

        ignore_time_budget = self.__get_debug_option('ignoreTimeBudget', False)
        if not ignore_time_budget:
            for budget_name, options in self.__time_budgets.items():
                budget_target_classes = set(options.target_classes)
                budget_variables = []
                # Collect all visits of targets that fall into to budget
                for (tidx, ivis), v in self.__variables.Tv_i.items():
                    target_class = self.__target_cache.target_class[tidx]
                    if target_class in budget_target_classes:
                        budget_variables.append(v)

                # TODO: This assumes a single, per visit exposure time which is fine for now
                #       but the logic could be extended further
                name = self.__make_name("Tv_i_sum", budget_name)
                constr = self.__visit_exp_time * self.__problem.sum([ v for v in budget_variables ]) <= 3600 * options.budget
                self.__constraints.Tv_i_sum[budget_name] = constr
                self.__add_constraint(name, constr)

    def __create_cobra_group_constraints(self, nvisits):
        # Make sure that there are enough targets in every cobra group for each visit

        ignore_cobra_group_minimum = self.__get_debug_option('ignoreCobraGroupMinimum', False)
        ignore_cobra_group_maximum = self.__get_debug_option('ignoreCobraGroupMaximum', False)

        for name, options in self.__cobra_groups.items():
            need_min = not ignore_cobra_group_minimum and options.min_targets is not None
            need_max = not ignore_cobra_group_maximum and options.max_targets is not None
            if need_min or need_max:                
                for ivis in range(nvisits):
                    for gidx in range(options.ngroups):
                        variables = self.__variables.CG_i[(name, ivis, gidx)]
                        if len(variables) > 0:
                            if need_min:
                                name = self.__make_name("Cv_CG_min", name, ivis, gidx)
                                constr = self.__problem.sum([ v for v in variables ]) >= options.min_targets
                                self.__constraints.CG_min[name, ivis, gidx] = constr
                                self.__add_constraint(name, constr)

                            if need_max:
                                name = self.__make_name("Cv_CG_max", name, ivis, gidx)
                                constr = self.__problem.sum([ v for v in variables ]) <= options.max_targets
                                self.__constraints.CG_max[name, ivis, gidx] = constr
                                self.__add_constraint(name, constr)

    def __create_unassigned_fiber_constraints(self, nvisits):
        # Make sure that enough fibers are kept unassigned, if this was requested
        # This is done by setting an upper limit on the sum of Cv_i edges

        ignore_reserved_fibers = self.__get_debug_option('ignoreReservedFibers', False)
        num_reserved_fibers = self.__get_netflow_option('numReservedFibers', 0)

        if not ignore_reserved_fibers and num_reserved_fibers > 0:
            max_assigned_fibers = self.__bench.cobras.nCobras - num_reserved_fibers
            for ivis in range(nvisits):
                variables = [var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == ivis]
                variables = [item for sublist in variables for item in sublist]

                name = self.__make_name("Cv_i_max", ivis)
                constr = self.__problem.sum(variables) <= max_assigned_fibers
                self.__constraints.Cv_i_max[ivis] = constr
                self.__add_constraint(name, constr)

    #endregion

    def solve(self):
        self.__problem.solve()
        self.__extract_assignments()

    def __extract_assignments(self):
        """Extract the fiber assignments from an LP solution"""

        nvisits = len(self.__visits)
        ncobras = self.__bench.cobras.nCobras

        self.__target_assignments = [ {} for _ in range(nvisits) ]
        self.__cobra_assignments = [ np.full(ncobras, -1, dtype=int) for _ in range(nvisits) ]

        def set_assigment(tidx, cidx, ivis):
            self.__target_assignments[ivis][tidx] = cidx
            self.__cobra_assignments[ivis][cidx] = tidx

        if self.__variables is not None:
            # This only works when the netflow problem is built
            for (tidx, ivis), vars in self.__variables.Tv_o.items():
                for (f, cidx) in vars:
                    if self.__problem.get_value(f) > 0:
                        set_assigment(tidx, cidx, ivis)
        else:
            # This also works when only the LP problem is loaded
            for k1, v1 in self.__problem._variables.items():
              if k1.startswith("Tv_Cv_"):
                    if self.__problem.get_value(v1) > 0:
                        _, _, tidx, cidx, ivis = k1.split("_")
                        set_assigment(int(tidx), int(cidx), int(ivis))

    def __extract_missed_science_targets(self):
        """Find science targets for each science target class that have been missed."""

        self.__missed_targets = { k: [] for k, options in self.__target_classes.items() if options.prefix == 'sci' }

        if self.__variables is not None:
            for tidx, f in self.__variables.T_i.items():
                if self.__problem.get_value(f) == 0:
                    target_class = self.__target_cache.target_class[tidx]
                    self.__missed_targets[target_class].append(tidx)
        else:
            raise NotImplementedError
        
    def __extract_partially_observed_science_targets(self):
        """Find science targets that are allocated to fibers during some visits but the total
        number of visits doesn't meet the requirements."""

        self.__partially_observed_targets = { k: [] for k, options in self.__target_classes.items() if options.prefix == 'sci' }

        if self.__variables is not None:
            for tidx, f in self.__variables.T_sink.items():
                if self.__problem.get_value(f) > 0:
                    target_class = self.__target_cache.target_class[tidx]
                    self.__partially_observed_targets[target_class].append(tidx)
        else:
            raise NotImplementedError()

    def __get_idcol(self, catalog):
        for c in ['objid', 'skyid', 'calid']:
            if c in catalog.data.columns:
                return c
        
        raise RuntimeError()
    
    def get_cobra_assignments(self):
        """
        Return a list of arrays for each visit that contain an array with the the associated target
        index to each cobra.
        """

        return self.__cobra_assignments
    
    def get_target_assignments(self, catalog):
        """
        Return a list of data frames for each visit with two columns: the object IDs and the fiber ID.
        """

        idcol = self.__get_idcol(catalog)

        assignments = []
        for i, p in enumerate(self.__visits):
            idx = np.array([ k for k in self.__target_assignments[i] ])                   # Target indices

            targets = pd.DataFrame(
                    {
                        f'target_{idcol}': [ np.int64(self.__targets.iloc[ti]['id']) for ti in idx ],
                        'fiberid': np.array([ self.__target_assignments[i][k] for k in idx ]) 
                    }
                ).astype({ 'fiberid': pd.Int32Dtype() }).set_index(f'target_{idcol}')

            # Find the index of each object in the catalog and return the matching fiberid
            assign = catalog.data[[idcol]].set_index(idcol).join(targets, how='left').reset_index()
            assignments.append(assign)

        return assignments

    def get_target_assignments_masks(self, catalog):
        """Returns a mask that indexes the selected targets within a catalog."""

        idcol = self.__get_idcol(catalog)
        assignments = self.get_target_assignments(catalog)
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

    def get_problematic_science_targets(self):
        if self.__missed_targets is None:
            self.__extract_missed_science_targets()

        if self.__partially_observed_targets is None:
            self.__extract_partially_observed_science_targets()

        return self.__missed_targets, self.__partially_observed_targets