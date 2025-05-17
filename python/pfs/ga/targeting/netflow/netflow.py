import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import pickle
from types import SimpleNamespace
from scipy.spatial import KDTree

from pfs.datamodel import TargetType, FiberStatus
# from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.CollisionSimulator2 import CollisionSimulator2
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID
from ics.cobraCharmer.cobraCoach import engineer

from .setup_logger import logger

from ..util.args import *
from ..util.config import *
from ..util.pandas import *
from ..data import Catalog
from ..instrument import CobraAngleFlags
from .visit import Visit
from .gurobiproblem import GurobiProblem
from .netflowexception import NetflowException
from .collisionsimulator import CollisionSimulator

class Netflow():
    """
    Implements the Network Flow algorithm to optimize fiber allocation. The problem
    is translated into an integer linear program (ILP)  solved by Gurobi.

    The initialization options are the following:

    solver_options : dict
        Options for the solver, passed to the GurobiProblem (or different class) constructor.
        See the documentation of the the LP solver API for details.
    netflow_options : dict
        Configuration options for the Netflow algorithm. The following options are supported:
        * black_dot_penalty : callable
            A function that returns the penalty for a cobra being too close to a black dot,
            as a function of distance. None if no penalty for being close to a black dot.
        * fiber_non_allocation_cost : float
            Cost for not allocating a fiber to a science or calibration target, or sky.
        * collision_distance : float
            Minimum distance between two fibers (ends or elbows) to avoid collision.
        * forbidden_targets : list
            List of forbidden target IDs. These are added to the problem as constraints.
        * forbidden_pairs : list of list
            List of forbidden target id pairs, where each pair is a list of two target IDs.
            These are added to the problem as constraints.
        * target_classes : dict
            Dictionary of target classes, see below for details.
        * cobra_groups : dict
            Dictionary of cobra groups, see below for details.
        * cobra_move_cost : callable
            A function that returns the cost of moving a cobra from the center,
            as a function of the distance. None if there is no extra cost.
        * time_budgets : dict
            Dictionary of time budgets for each target class, see below for details.
        * num_reserved_fibers : int
            Number of reserved fibers, without any additional constraints.
        * allow_more_visits : bool
            If True, allow more than visits per science targets than required. Default is False.
        * epoch : float
            Epoch for proper motion and parallax calculations, all catalogs must be at the same epoch.
        * ignore_proper_motion : bool
            If True, ignore proper motion when calculating the focal plane positions.
        * fiberids_path : str
            Path to the fiberids file. Should point to `pfs_utils/data/fiberids`. 
    debug_options : dict
        Debugging options for the Netflow algorithm. When an option is set to True, the
        corresponding constraint are created but not added to the problem. The following
        options are supported:
        * ignore_endpoint_collisions
        * ignore_elbow_collisions
        * ignore_broken_cobra_collisions
        * ignore_forbidden_targets
        * ignore_forbidden_pairs
        * ignore_calib_target_class_minimum
        * ignore_calib_target_class_maximum
        * ignore_science_target_class_minimum
        * ignore_science_target_class_maximum
        * ignore_time_budget
        * ignore_cobra_group_minimum
        * ignore_cobra_group_maximum
        * ignore_reserved_fibers
    targetClasses : dict
        Dictionary of target classes. The keys are the names of the target classes. Two special
        keys are reserved: 'sky' and 'cal', for sky fibers and calibration targets, respectively.
        Each target class definition is a dictionary with the following keys
        * prefix : str
            Prefix of the target class, e.g. 'sci', 'cal', 'sky'. This is used to identify the
            type of the target class, whereas names can be arbitrary.
        * min_targets : int
            Minimum number of targets of this class that must be observed.
        * max_targets : int
            Maximum number of targets of this class that can be observed.
        * non_observation_cost : float
            Cost for not observing a target of this class.
        * partial_observation_cost : float
            Cost for partially observing a target of this class. Must be larger than
            the non-obervation cost.
    cobraGroups : dict
        Dictionary of cobra groups. The keys are the names of the cobra groups. Each cobra group
        definition is a dictionary with the following keys:
        * groups : ndarray
            An array of integers indicating the group of each cobra.
        * target_classes : list
            A list of target classes that should be considered when creating the constraints
            applicable to this cobra group.
        * min_targets : int
            Minimum number of targets of the target classes that must be observed with this cobra group.
        * max_targets : int
            Maximum number of targets of the target classes that can be observed with this cobra group.
        * non_observation_cost : float
            Cost for not observing a target of any of the target classes with this cobra group.

    After initialization, the targets are added in form of `Catalog` objects which in
    turn contain a Pandas dataframe in the variables `data`. The targets are added
    with the function `append_targets` and its wrappers. The input dataframe is converted
    into a uniform format with the following columns:
    * targetid     - must be unique across all targets
    * RA          - right ascension in degrees
    * Dec         - declination in degrees
    * penalty     - penalty for observing the target, smaller for higher priory targets
    * prefix      - target prefix, e.g. 'sci', 'cal', 'sky'
    * exp_time    - total requested exposure time
    * priority    - priority of the target, only used to generate a name for the target class e.g. 'sci_P1'.
    * class       - name of the target class to which the target belongs
    
    Using dataframes helps with the manipulation of the data and prevents instantiation a
    large number of objects.
    
    Once the targets are added, the problem is built with the function `build`. This function
    constructs the ILP problem by defining the variables and constraints. First, the number of
    required visits for each target it calculated in `__calculate_target_visits`. The duration
    of each visit is assumed to be the same for all pointings and visits. Next, the targets are
    cached, which consists of rewriting the pandas dataframe `__targets` into numpy arrays by
    the function `__cache_targets`. The advantage of this is the much faster item access by index
    compared to the pandas dataframe. Everywhere where the targets are indexed with variables
    named similar to `target_idx` or `tidx`, they're indices into the arrays of  `__target_cache`.
    Note that `target_idx` is different from `target_id`, where the later is the unique identifier
    of the target within the catalog, as provided by the user.
    
    A `Visit` object is created for each pointing and visit. These are stored in the list `__visits`.
    Everywhere where variables named similar to `visit_idx` or `vidx` are used, they're indices into
    the list `__visits`.
    
    In the next step, the focal plane coordinates of the targets are calculated for each pointing by
    `__calculate_target_fp_pos`. Visits of the same pointing will result in the same focal plane positions.
    
    After the setup, the ILP problem is built in `__build_ilp_problem`. This is the second most time-
    consuming step of the algorithm, after the actual optimization, since it requires creating
    millions of variables and constraints. After the ILP problem is built it can be saved to a file
    with `save_problem` and loaded with `load_problem`.
    
    The problem can be solved with the `solve` function which executes the solver and
    `extract_assignments` can be called to extract the results. The results are stored in the variables
    `__target_assignments` and `__cobra_assignments` which are lists of dictionaries keyed by
    the target index or the cobra index, respectively.

    Once the problem is solved and target assignments 
    """
    
    def __init__(self, 
                 name,
                 instrument,
                 pointings,
                 workdir=None,
                 filename_prefix=None,
                 solver_options=None,
                 netflow_options=None,
                 debug_options=None):

        # Configurable options

        self.__name = name                          # Problem name
        self.__instrument = instrument              # Instrument object
        self.__pointings = pointings                # List of pointings
        self.__visits = None                        # List of visits
        
        self.__workdir = workdir if workdir is not None else os.getcwd()
        self.__filename_prefix = filename_prefix if filename_prefix is not None else ''
        self.__solver_options = solver_options
        self.__netflow_options = netflow_options
        self.__debug_options = debug_options
        self.__resume = False                        # Resume flag

        # Internally used variables

        self.__use_named_variables = None

        self.__bench = instrument.bench             # Bench object
        self.__fiber_map = instrument.fiber_map     # Grand fiber map
        self.__blocked_fibers = instrument.blocked_fibers
        self.__calib_model = instrument.calib_model

        self.__problem_type = GurobiProblem
        self.__problem = None                       # ILP problem, already wrapped

        self.__name_counter = None                  # Variable name counter

        self.__variables = None                     # LP variables, see __init_variables
        self.__constraints = None                   # LP constraints, see __init_constraints

        self.__forbidden_targets = None             # List of forbidden individual target, identified by target_idx
        self.__forbidden_pairs = None               # List of forbidden target pairs, identified by target_idx
        self.__black_dots = None                    # Nearest black dots around each cobra

        self.__visit_exp_time = None                # Exposure time of a single visit in integer seconds
        self.__targets = None                       # DataFrame of targets
        self.__target_cache = None                  # Cached targets in numpy arrays
        self.__target_fp_pos = None                 # Target focal plane positions for each pointing

        self.__target_mask = None                   # Mask for valid targets for each pointing
        self.__target_cache_mask = None             # Mask for targets in the target cache
        self.__target_to_cache_map = None           # Maps target index to cache index
        self.__cache_to_target_map = None           # Maps cache index to target index
        self.__cache_to_fp_pos_map = None           # Maps tidx to fpidx for each pointing
        self.__fp_pos_to_cache_map = None           # Maps fpidx to tidx for each pointing

        self.__visibility = None                    # Visibility and elbow positions for each target for each pointing
        self.__collisions = None                    # Tip and elbow collision matrices for each pointing
        
        self.__target_assignments = None            # List of dicts for each visit, keyed by target index
        self.__cobra_assignments = None             # List of dicts for each visit, keyed by cobra index

        self.__cobra_collisions = None              # List of arrays for each visit, keyed by cobra index

    #region Property accessors

    def __get_name(self):
        return self.__name
    
    def __set_name(self, value):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_netflow_options(self):
        return self.__netflow_options

    netflow_options = property(__get_netflow_options)

    def __get_solver_options(self):
        return self.__solver_options
    
    solver_options = property(__get_solver_options)

    def __debug_options(self):
        return self.__debug_options
    
    debug_options = property(__debug_options)

    def __get_pointings(self):
        return self.__pointings
    
    pointings = property(__get_pointings)

    def __get_visits(self):
        return self.__visits
    
    visits = property(__get_visits)

    def __get_targets(self):
        return self.__targets
    
    def _set_targets(self, value):
        self.__targets = value
    
    targets = property(__get_targets)

    def __get_target_classes(self):
        return self.__netflow_options.target_classes
    
    target_classes = property(__get_target_classes)

    def __get_target_assignments(self):
        return self.__target_assignments
    
    target_assignments = property(__get_target_assignments)

    def __get_cobra_assignments(self):
        return self.__cobra_assignments
    
    cobra_assignments = property(__get_cobra_assignments)

    def __get_cobra_collisions(self):
        return self.__cobra_collisions
    
    cobra_collisions = property(__get_cobra_collisions)

    #endregion

    def __validate_options(self):

        for name, options in self.__netflow_options.target_classes.items():                
            # Sanity check for science targets: make sure that partialObservationCost
            # is larger than nonObservationCost
            if options.prefix in self.__netflow_options.science_prefix \
                and options.non_observation_cost is not None \
                and options.partial_observation_cost is not None \
                and options.partial_observation_cost < options.non_observation_cost:
                
                raise NetflowException(
                    f"Found target class `{name}` where partial_observation_cost "
                    "is smaller than non_observation_cost")

    def __get_netflow_option(self, value, default=None):
        """
        Substitute the default value if value is None. This function
        is kept for compatibility after switching to config classes
        from simple dictionaries.
        """

        if value is not None:
            return value
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
    
    def __get_targets_filename(self, filename=None):
        # TODO: add support for .zip, .gz, .bz2, and .7z
        fn = self.__get_prefixed_filename(filename, 'netflow_targets.feather')
        fn = self.__append_extension(fn, '.feather')
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
        """Save the targets in a file"""
        
        fn = self.__get_targets_filename(filename)
        self.__targets.to_feather(fn)

        logger.debug(f'Saved targets to `{fn}`.')

    def load_targets(self, filename=None):
        """Load the targets from a file"""

        fn = self.__get_targets_filename(filename)
        self.__targets = pd.read_feather(fn)

        logger.debug(f'Loaded targets from `{fn}`.')

    def load_problem(self, filename=None, options=None):
        """Read the LP problem from a file"""

        fn = self.__get_problem_filename(filename)
        if os.path.exists(fn):
            self.__problem = self.__problem_type(self.__name, self.__solver_options)
            self.__problem.read_problem(fn)
        else:
            raise FileNotFoundError(f"Problem file `{fn}` not found.")

    def __restore_variables(self):

        logger.info("Restoring variables...")

        self.__init_variables()        
        vars = self.__problem.get_variables()

        for i, f in enumerate(vars):
            if f.varName == 'cost':
                continue
            elif f.varName.startswith('STC_sink'):        # STC_sink_sci_P0
                m = re.match(r'STC_sink_(\w+)', f.varName)
                target_class = m.group(1)
                self.__variables.STC_sink[target_class] = f
            elif f.varName.startswith('STC_T'):           # STC_T[202700]
                m = re.match(r'STC_T\[(\d+)\]', f.varName)
                tidx = int(m.group(1))
                target_class = self.__target_cache.target_class[tidx]
                self.__variables.T_i[tidx] = f
                self.__variables.STC_o[target_class].append(f)
            elif f.varName.startswith('CTCv_sink'):     # CTCv_sink_sky_0
                m = re.match(r'CTCv_sink_(\w+)_(\d+)', f.varName)
                target_class = m.group(1)
                vidx = int(m.group(2))
                self.__variables.CTCv_sink[(target_class, vidx)] = f
            elif f.varName.startswith('CTCv_Tv'):       # CTCv_Tv_0[22810]
                m = re.match(r'CTCv_Tv_(\d+)\[(\d+)\]', f.varName)
                vidx = int(m.group(1))
                tidx = int(m.group(2))
                target_class = self.__target_cache.target_class[tidx]
                self.__variables.Tv_i[(tidx, vidx)] = f
                self.__variables.CTCv_o[(target_class, vidx)].append(f)
            elif f.varName.startswith('T_sink'):        # T_sink[202611]
                m = re.match(r'T_sink\[(\d+)\]', f.varName)
                tidx = int(m.group(1))
                self.__variables.T_sink[tidx] = f
            elif f.varName.startswith('T_Tv'):          # T_Tv_0[249940]
                m = re.match(r'T_Tv_(\d+)\[(\d+)\]', f.varName)
                vidx = int(m.group(1))
                tidx = int(m.group(2))
                self.__variables.T_o[tidx].append(f)
                self.__variables.Tv_i[(tidx, vidx)] = f
            elif f.varName.startswith('Tv_Cv'):         # Tv_Cv_0_0[22810]
                m = re.match(r'Tv_Cv_(\d+)_(\d+)\[(\d+)\]', f.varName)
                vidx = int(m.group(1))
                cidx = int(m.group(2))
                tidx = int(m.group(3))
                target_class = self.__target_cache.target_class[tidx]
                self.__variables.Cv_i[(cidx, vidx)].append((f, tidx))
                self.__variables.Tv_o[(tidx, vidx)].append((f, cidx))
                for cg_name, options in self.__netflow_options.cobra_groups.items():
                    if target_class in options.target_classes:
                        self.__variables.CG_i[(cg_name, vidx, options.groups[cidx])].append(f)
            else:
                raise NotImplementedError()
            
            self.__variables.all[f.varName] = f

        logger.info(f'Restored {len(vars)} variables.')

    def __restore_constraints(self):

        logger.info("Restoring constraints...")

        self.__init_constraints()        
        constrs = self.__problem.get_constraints()

        for i, c in enumerate(constrs):
            if c.constrName.startswith('Tv_o_coll'):
                m = re.match(r'^Tv_o_coll_(\d+)_(\d+)_(\d+)$', c.constrName)
                if m is not None:
                    tidx1 = int(m.group(1))
                    tidx2 = int(m.group(2))
                    vidx = int(m.group(3))
                    self.__constraints.Tv_o_coll[(tidx1, tidx2, vidx)] = c
                else:
                    m = re.match(r'^Tv_o_coll_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)$', c.constrName)
                    tidx1 = int(m.group(1))
                    cidx1 = int(m.group(2))
                    tidx2 = int(m.group(3))
                    cidx2 = int(m.group(4))
                    vidx = int(m.group(5))
                    self.__constraints.Tv_o_coll[(tidx1, cidx1, tidx2, cidx2, vidx)] = c
            elif c.constrName.startswith('Tv_o_broken'):
                m = re.match(r'^Tv_o_broken_(\d+)_(\d+)_(\d+)_(\d+)$', c.constrName)
                cidx = int(m.group(1))
                ncidx = int(m.group(2))
                tidx = int(m.group(3))
                vidx = int(m.group(4))
                self.__constraints.Tv_o_broken[(cidx, ncidx, tidx, vidx)] = c
            elif c.constrName.startswith('Tv_o_forb'):
                m = re.match(r'^Tv_o_forb_(\d+)_(\d+)_(\d+)$', c.constrName)
                tidx1 = int(m.group(1))
                tidx2 = int(m.group(2))
                vidx = int(m.group(3))
                self.__constraints.Tv_o_forb[(tidx1, tidx2, vidx)] = c
            elif c.constrName.startswith('STC_o_sum'):
                m = re.match(r'^STC_o_sum_(\w+)$', c.constrName)
                target_class = m.group(1)
                self.__constraints.STC_o_sum[target_class] = c
            elif c.constrName.startswith('STC_o_min'):
                m = re.match(r'^STC_o_min_(\w+)$', c.constrName)
                target_class = m.group(1)
                self.__constraints.STC_o_min[target_class] = c
            elif c.constrName.startswith('STC_o_max'):
                m = re.match(r'^STC_o_max_(\w+)$', c.constrName)
                target_class = m.group(1)
                self.__constraints.STC_o_max[target_class] = c
            elif c.constrName.startswith('CTCv_o_sum'):
                m = re.match(r'^CTCv_o_sum_(\w+)_(\d+)$', c.constrName)
                target_class = m.group(1)
                vidx = int(m.group(2))
                self.__constraints.CTCv_o_sum[(target_class, vidx)] = c
            elif c.constrName.startswith('CTCv_o_min'):
                m = re.match(r'^CTCv_o_min_(\w+)_(\d+)$', c.constrName)
                target_class = m.group(1)
                vidx = int(m.group(2))
                self.__constraints.CTCv_o_min[(target_class, vidx)] = c
            elif c.constrName.startswith('CTCv_o_max'):
                m = re.match(r'^CTCv_o_max_(\w+)_(\d+)$', c.constrName)
                target_class = m.group(1)
                vidx = int(m.group(2))
                self.__constraints.CTCv_o_max[(target_class, vidx)] = c
            elif c.constrName.startswith('T_i_T_o_sum_0'):
                m = re.match(r'^T_i_T_o_sum_0_(\d+)$', c.constrName)
                tidx = int(m.group(1))
                self.__constraints.T_i_T_o_sum[(tidx, 0)] = c
            elif c.constrName.startswith('T_i_T_o_sum_1'):
                m = re.match(r'^T_i_T_o_sum_1_(\d+)$', c.constrName)
                tidx = int(m.group(1))
                self.__constraints.T_i_T_o_sum[(tidx, 1)] = c
            elif c.constrName.startswith('T_i_T_o_sum'):
                m = re.match(r'^T_i_T_o_sum_(\d+)$', c.constrName)
                tidx = int(m.group(1))
                self.__constraints.T_i_T_o_sum[tidx] = c
            elif c.constrName.startswith('Tv_i_Tv_o_sum'):
                m = re.match(r'^Tv_i_Tv_o_sum_(\d+)_(\d+)$', c.constrName)
                tidx = int(m.group(1))
                vidx = int(m.group(2))
                self.__constraints.Tv_i_Tv_o_sum[(tidx, vidx)] = c
            elif c.constrName.startswith('Cv_i_sum'):
                m = re.match(r'^Cv_i_sum_(\d+)_(\d+)$', c.constrName)
                cidx = int(m.group(1))
                vidx = int(m.group(2))
                self.__constraints.Cv_i_sum[(cidx, vidx)] = c
            elif c.constrName.startswith('Tv_i_sum'):
                m = re.match(r'^Tv_i_sum_(\w+)$', c.constrName)
                budget_name = m.group(1)
                self.__constraints.Tv_i_sum[budget_name] = c
            elif c.constrName.startswith('Cv_CG_min'):
                m = re.match(r'^Cv_CG_min_(\w+)_(\d+)_(\d+)$', c.constrName)
                cg_name = m.group(1)
                vidx = int(m.group(2))
                cidx = int(m.group(3))
                self.__constraints.Cv_CG_min[(cg_name, vidx, cidx)] = c
            elif c.constrName.startswith('Cv_CG_max'):
                m = re.match(r'^Cv_CG_max_(\w+)_(\d+)_(\d+)$', c.constrName)
                cg_name = m.group(1)
                vidx = int(m.group(2))
                cidx = int(m.group(3))
                self.__constraints.Cv_CG_min[(cg_name, vidx, cidx)] = c
            elif c.constrName.startswith('Cv_i_max'):
                m = re.match(r'^Cv_i_max_(\d+)$', c.constrName)
                vidx = int(m.group(1))
                self.__constraints.Cv_i_max[vidx] = c
            else:
                raise NotImplementedError()

            self.__constraints.all[c.constrName] = c

        logger.info(f'Restored {len(constrs)} constraints.')

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
    
    def __get_forbidden_targets_config(self):
        """
        Look up the target index based on the target id of forbidden targets.
        """

        forbidden_targets = self.__netflow_options.forbidden_targets
        fpp = []
        wrong_id_count = 0
        if forbidden_targets is not None:
            for i, target_id in enumerate(forbidden_targets):            
                tidx = self.__target_to_cache_map[self.__targets.index.get_loc(target_id)]
                if tidx != -1:
                    fpp.append(tidx)
                else:
                    wrong_id_count += 1
                    
        if wrong_id_count > 0:
            logger.warning(f"Found {wrong_id_count} forbidden targets that are outside the pointings.")
                
        return fpp
    
    def __get_forbidden_pairs_config(self):
        """
        Look up the target indices based on the target ids of forbidden pairs.
        """

        forbidden_pairs = self.__netflow_options.forbidden_pairs
        fpp = []
        wrong_id_count = 0
        if forbidden_pairs is not None:
            for i, pair in enumerate(forbidden_pairs):
                if len(pair) != 2:
                    raise NetflowException(f"Found an incorrect number of target ids in forbidden pair list at index {i}.")
            
                tidx_list = [ self.__target_to_cache_map[self.__targets.index.get_loc(p)] for p in pair ]
                if -1 not in tidx_list:
                    fpp.append(tidx_list)
                else:
                    wrong_id_count += 1

        if wrong_id_count:
            logger.warning(f"Found {wrong_id_count} forbidden targets pairs that are outside the pointings.")
                
        return fpp

    def __get_black_dot_config(self):
        black_dot_penalty = self.__get_netflow_option(self.__netflow_options.black_dot_penalty, None)
        
        center, radius = self.__get_closest_black_dots()
        black_dots = SimpleNamespace(
            # Closest black dots for each cobra
            black_dot_center = center,
            black_dot_radius = radius,
            black_dot_penalty = black_dot_penalty
        )

        logger.info('Successfully loaded black dot configuration.')

        return black_dots
    
    #endregion
    #region PFI functions

    def __calculate_fp_pos(self, pointing,
                           ra, dec,
                           pmra=None, pmdec=None, parallax=None, rv=None, epoch=None,
                           ignore_proper_motion=False):

        if ignore_proper_motion:
            parallax = None
            pmra = pmdec = None
            rv = None

        fp_pos = self.__instrument.radec_to_fp_pos(ra, dec,
                                                   pointing=pointing,
                                                   epoch=epoch,
                                                   pmra=pmra, pmdec=pmdec, parallax=parallax, rv=rv)
                
        return fp_pos
    
    def __get_closest_black_dots(self):
        """
        For each cobra, return the list of focal plane black dot positions.
        
        Returns
        -------
        list:
            List of arrays of complex focal plane positions that include the black dots
            of the cobra and its neighbors.
        """

        # NOTE: the black dot radius multiplier instrument option `black_dot_radius_margin` is passed
        #       to the `Bench` object which already multiplies the radii by it so we do not need
        #       to multiple here again

        center = []
        radius = []
        for cidx in range(len(self.__bench.cobras.centers)):
            nb = self.__bench.getCobraNeighbors(cidx)
            center.append(self.__bench.blackDots.centers[[cidx] + list(nb)])
            radius.append(self.__bench.blackDots.radius[[cidx] + list(nb)])

        return center, radius

    def __get_closest_black_dot_distance(self, cidx, fp_pos):
        """
        Calculate the distance of a target from the closest black dot. Only the
        black dot closest to the cobra and its neighbors are considered, as in
        `__get_closest_black_dots`
        """

        fp_pos = np.atleast_1d(fp_pos)

        # Distance from each nearby black dot
        dist = np.abs(self.__black_dots.black_dot_center[cidx] - fp_pos[:, None])
        
        # Take minimum distance from each target
        return np.min(dist, axis=-1)

    def __get_cobra_center_distance(self, cidx, fp_pos):
        """
        Calculate the distance of the target from the center of the cobra, in focal plane
        coordinates.
        """

        fp_pos = np.atleast_1d(fp_pos)
        return np.abs(self.__bench.cobras.centers[cidx] - fp_pos)
        
    def __build_fp_pos_kdtree(self, fp_pos):
        # Build a KD-tree on the target focal plane coordinates

        leaf_size = 3 * fp_pos.size // self.__bench.cobras.nCobras
        leaf_size = max(2, leaf_size)

        # Construct the KD tree
        kdtree = KDTree(np.column_stack((fp_pos.real, fp_pos.imag)), leafsize=leaf_size)
        
        return kdtree
   
    def __calculate_visibility_and_elbow(self, fp_pos, fpidx_to_tidx_map):
        """
        Calculate the visibility and the corresponding elbow position for each
        cobra from the focal plane positions in ˙fp_pos˙.

        This is similar to what's implemented in ics.cobraOps.TargetSelector but
        supports better vectorization. It also calculates the cobra motor angles and
        elbow positions and filters out invalid solutions.
        """

        # Create a KDTree to find targets within the patrol radisu
        kdtree = self.__build_fp_pos_kdtree(fp_pos)

        # Query the tree for each cobra in parallel
        centers = self.__instrument.bench.cobras.centers
        rmin = self.__instrument.bench.cobras.rMin
        rmax = self.__instrument.bench.cobras.rMax
        
        # Do a coarse search first because we'll have to calculate the distances anyway
        idx = kdtree.query_ball_point(np.stack([centers.real, centers.imag], axis=1), rmax, eps=rmax.max() * 0.05)

        # Build two dictionaries with different look-up directions
        targets_cobras = defaultdict(list)
        cobras_targets = defaultdict(list)

        for cidx in range(self.__bench.cobras.nCobras):
            # Skip broken cobra
            if self.__bench.cobras.hasProblem[cidx]:
                continue

            # Calculate the cobra angles and elbow positions
            fpidx = np.array(idx[cidx], dtype=int)
            if fpidx.size > 0:
                # Filter out targets too close to any of the black dots
                if self.__black_dots is not None:
                    black_dot_center = self.__black_dots.black_dot_center[cidx]
                    black_dot_radius = self.__black_dots.black_dot_radius[cidx]
                    black_dot_dist = np.sqrt(np.abs(fp_pos[fpidx] - black_dot_center[..., None]))

                    # Target must be far from all black nearby black dots
                    black_dot_mask = np.all(black_dot_dist > black_dot_radius[..., None], axis=0)
                else:
                    black_dot_mask = None

                # Find the angles and elbow positions for the given focal plane coordinates
                theta, phi, d, eb_pos, flags = self.__instrument.fp_pos_to_cobra_angles(fp_pos[fpidx], np.atleast_1d(cidx))

                # There can be up to two solutions for the cobra angles, given a single
                # focal plane position. Filter out invalid solutions and collect elbow positions
                elbows = defaultdict(tuple)
                angles = defaultdict(tuple)
                for i in range(2):
                    mask = (flags[..., i] & CobraAngleFlags.SOLUTION_OK) != 0
                    mask &= (flags[..., i] & (CobraAngleFlags.TOO_CLOSE_TO_CENTER | CobraAngleFlags.TOO_FAR_FROM_CENTER)) == 0
                    mask &= (flags[..., i] & (CobraAngleFlags.THETA_OUT_OF_RANGE | CobraAngleFlags.PHI_OUT_OF_RANGE)) == 0

                    if black_dot_mask is not None:
                        mask &= black_dot_mask

                    # Map focal plane index to the target cache index
                    tidx = fpidx_to_tidx_map[fpidx[mask]]
                    ebp = eb_pos[..., i][mask]
                    th = theta[..., i][mask]
                    ph = phi[..., i][mask]

                    for ti in range(tidx.size):
                        elbows[tidx[ti]] = elbows[tidx[ti]] + (ebp[ti],)
                        angles[tidx[ti]] = angles[tidx[ti]] + ((th[ti], ph[ti]),)
                    
                for ti in elbows:
                    targets_cobras[ti].append((cidx, elbows[ti], angles[ti]))
                    cobras_targets[cidx].append((ti, elbows[ti], angles[ti]))

        return targets_cobras, cobras_targets

    def __get_targets_in_intersections(self, pidx, cidx):

        visibility = self.__visibility[pidx]

        # Determine target indices visible by this cobra and its neighbors
        tidx = set(ti for ti, _, _ in visibility.cobras_targets[cidx])

        # Neighboring cobras
        ntidx = {}
        for ncidx in self.__bench.getCobraNeighbors(cidx):
            if ncidx in visibility.cobras_targets:
                ntidx[ncidx] = np.array(list(set(ti for ti, _, _ in visibility.cobras_targets[ncidx]).intersection(tidx)), dtype=int)

        return ntidx
    
    def __get_colliding_endpoints(self, pidx, collision_distance):
        """Return the list of target pairs that would cause fiber
        collision when observed by two neighboring fibers.
        
        Parameters
        ==========

        pidx : int
            Index of the pointing
        collision_distance:
            Maximum distance causing a collision.

        Returns
        =======
        set : pairs of target indices that would cause fiber tip collisions.
        """

        fp_pos = self.__target_fp_pos[pidx]
        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]
        visibility = self.__visibility[pidx]

        pairs = set()
        for cidx in visibility.cobras_targets:
            # Determine target indices visible by this cobra and its neighbors
            intersections = self.__get_targets_in_intersections(pidx, cidx)

            for ncidx, ntidx in intersections.items():
                # Calculate the distance matrix of focal plane positions
                fp = fp_pos[tidx_to_fpidx_map[ntidx]]
                d = np.abs(np.subtract.outer(fp, fp))

                # Generate the list of pairs
                for m, n in zip(*np.where(d < collision_distance)):
                    # Only store pairs once
                    if m < n:
                        pairs.add((ntidx[m], ntidx[n]))

        return pairs
    
    def __get_colliding_elbows(self, pidx, collision_distance):
        """
        For each target-cobra pair, and the corresponding elbow position,
        return the list of other targets that are too close to the "upper arm" of
        the cobra.
        
        Parameters
        ==========

        pidx : int
            Index of the pointing
        collision_distance:
            Maximum distance causing a collision.

        Returns
        =======
        dict : Dictionary of list of targets too close to the cobra indexed by
               all possible target-cobra pairs.
        """

        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]
        visibility = self.__visibility[pidx]
        fp_pos = self.__target_fp_pos[pidx]

        res = defaultdict(list)
        for cidx, tidx_elbows_angles in visibility.cobras_targets.items():            
            # Determine target indices visible by neighbors of this cobra
            for ncidx in self.__bench.getCobraNeighbors(cidx):
                if ncidx in visibility.cobras_targets:
                    ntidx_elbows_angles = visibility.cobras_targets[ncidx]

                    # Loop over targets visible by the first cobra
                    for tidx, elbows, angles in tidx_elbows_angles:
                        fp = fp_pos[[tidx_to_fpidx_map[tidx]]]

                        # Check all possible elbow positions but make sure one collision is counted
                        # only once, otherwise we get duplicate (tidx1, cidx1, tidx2, cidx2, vidx) values
                        # TODO: now we exclude a target if ANY of the angle solution collide with a neighbor
                        #       but this could be changed to an OR
                        ntidx = []
                        for eb in elbows:
                            ntidx += [ ti for ti, _, _ in ntidx_elbows_angles ]
                        
                        ntidx = np.unique(ntidx)
                        nfp = fp_pos[tidx_to_fpidx_map[ntidx]]
                        d = self.__bench.distancesToLineSegments(nfp, np.repeat(fp, ntidx.shape), np.repeat(eb, ntidx.shape))
                        for nti in ntidx[d < collision_distance]:
                            res[(cidx, tidx)].append((ncidx, nti))

        return res
    
    def __get_colliding_broken_cobras(self, pidx, collision_distance):
        """
        For each broken cobra, look up the neighboring cobras and return the list of targets
        that would cause endpoint or elbow collision with the broken cobra in home position.

        The collision detection is done the same way as in CollisionSimulator and CollisionSimulator2.
        """

        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]
        visibility = self.__visibility[pidx]

        # For each broken cobra, look up the neighboring cobras and loop over them
        res = defaultdict(list)
        for cidx in np.arange(self.__bench.cobras.nCobras)[self.__bench.cobras.hasProblem]:
            fp_pos = np.atleast_1d(self.__bench.cobras.home0[cidx])
            eb_pos = self.__bench.cobras.calculateCobraElbowPositions(cidx, fp_pos)

            for ncidx in self.__bench.getCobraNeighbors(cidx):
                # This is another broken cobra nothing to do with it
                if self.__bench.cobras.hasProblem[ncidx]:
                    continue
                
                # Find the targets that are visible by the neighboring cobra
                ntidx = np.array([ ti for ti, _, _ in visibility.cobras_targets[ncidx] ])
                
                if ntidx.size > 0:
                    nfp_pos = self.__target_fp_pos[pidx][tidx_to_fpidx_map[ntidx]]

                    # TODO: elbows and angles might have multiple elements if there are multiple solutions
                    #       for the cobra angles. Here we only use the first solution
                    neb_pos = np.array([ eb[0] for ti, eb, an in visibility.cobras_targets[ncidx] ])

                    # Distance from the home position of the broken cobra to the targets
                    d = self.__instrument.bench.distancesBetweenLineSegments(
                        np.repeat(fp_pos, ntidx.size),
                        np.repeat(eb_pos, ntidx.size),
                        nfp_pos,
                        neb_pos)

                    t = list(ntidx[d < collision_distance])
                    if len(t) > 0:
                        res[(cidx, ncidx)] += t

        return res
    
    def __get_colliding_trajectories(self, pidx, collision_distance):
        """
        Calculate the trajectories of cobras belonging to targets in the overlap regions.

        This is based on the class CollisionSimulator2 but allows for better vectorization
        """

        # TODO: bring out parameters
        time_step = 20
        max_steps = 2000

        # Initialize the cobra coach engineer
        engineer.setCobraCoach(self.__instrument.cobra_coach)
        engineer.setConstantOntimeMode(maxSteps=max_steps)

        # Collect all targets that are in the overlap regions
        # TODO: this could be cached at some point
        # TODO: overlap region should be revised because cobra tip and elbow can be further away than
        #       the fiber

        # TODO: we can exclude the colliding tips and elbows from the checks

        fp_pos = self.__target_fp_pos[pidx]
        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]

        def calcuate_trajectory(cidx, tidx):
            theta, phi, d, eb_pos, flags = self.__instrument.fp_pos_to_cobra_angles(
                fp_pos[tidx_to_fpidx_map[tidx]], np.full_like(tidx, cidx))
            
            raise NotImplementedError()

            return fp, eb
        
        for cidx in range(self.__bench.cobras.nCobras):
            intersections = self.__get_targets_in_intersections(pidx, cidx)

            # Calculate the trajectory for each focal plane position
            for ncidx, tidx in intersections.items():
                # Calculate the trajectory for the cobras
                fp1, eb1 = calcuate_trajectory(cidx, tidx)
                fp2, eb2 = calcuate_trajectory(ncidx, tidx)
                
        pass

    #endregion

    def __filter_targets(self, data: pd.DataFrame, pointings):
        # Filter down the data to the field of view, this saves a lot of time
        # when calculating the focal plane positions
        mask = []
        cc = SkyCoord(np.array(data['RA']) * u.deg, np.array(data['Dec']) * u.deg)
        for p in pointings:
            pp = SkyCoord(p.ra * u.deg, p.dec * u.deg)
            sep = cc.separation(pp).degree
            mask.append(sep < 0.75)

        return mask

    def __discover_id_column(self, prefix, df):
        # Try to discover the identity column in a data frame
        
        if 'targetid' in df:
            return 'targetid'
        elif 'objid' in df:
            return 'objid'
        elif 'skyid' in df:
            return 'skyid'
        elif prefix == 'sci':
            return 'objid'
        elif prefix == 'cal':  
            return 'objid'
        elif prefix == 'sky':
            return 'skyid'
        else:
            return None

    def get_catalog_data(self,
                         catalog,
                         prefix,
                         exp_time=None,
                         priority=None,
                         penalty=None,
                         non_observation_cost=None,
                         mask=None,
                         filter=None,
                         selection=None):
        
        """
        Get a dataframe containing the columns of an input target list. Some
        columns are updated or generated based on the input parameters and
        the netflow configuration.

        By default, use values such as `exp_time`, `priority` and `penalty` from
        the input catalog but when these values are passed to the function, they
        will override the values in the catalog.

        Parameters
        ----------
        catalog : Catalog or pd.DataFrame
            The input catalog to be added to the target list. If a DataFrame is passed,
            it must contain the columns `RA`, `Dec` and `__key`.
        prefix : str
            The prefix to be used for the target class. This is usually `sci` for science
            targets and `cal` for calibration targets.
        exp_time : float, optional
            The exposure time to be used for the targets. If not specified, the value
            from the input catalog will be used.
        priority : int, optional
            The priority to be used for the targets. If not specified, the value
            from the input catalog will be used.
        penalty : int, optional
            The penalty to be used for the targets. If not specified, the value
            from the input catalog will be used.
        non_observation_cost : int, optional
            The non-observation cost to be used for the targets. If not specified,
            the value from the input catalog will be used.
        mask : array-like, optional
            A mask to filter the input catalog. If not specified, all targets will
            be used.
        filter : callable, optional
            A callable that acts as filter to be applied to the input catalog.
            If not specified, no filter will be applied. See `Observation.get_data` for details.
        selection : Selection, optional
            A selection object to be applied to the input catalog. If not specified,
            no selection will be applied. See `Observation.get_data` for details.
        """

        # Depending on the type of catalog, we need to filter the data
        # differently. For a DataFrame, we can use the mask directly, but for a
        # Catalog, we need to use the `get_data` method to filter the data.
        if isinstance(catalog, Catalog):
            df = catalog.get_data(mask=mask, filter=filter, selection=selection)
        elif isinstance(catalog, pd.DataFrame):
            if filter is not None:
                raise NetflowException("The `filter` parameter is not supported for DataFrames.")
            
            if selection is not None:
                raise NetflowException("The `selection` parameter is not supported for DataFrames.")
            
            df = catalog[mask if mask is not None else ()]
        
        # Identify the identity column in the input catalog
        id_column = self.__discover_id_column(prefix, df)
        if id_column is None:
            raise NetflowException(f"Cannot determine identity column for catalog with prefix `{prefix}`.")   
        
        # Create the formatted dataset from the input catalog
        # Make sure that we convert to the right data type everywhere
        targets = pd.DataFrame({
            'key': df['__key'].astype(str),
            'targetid': df[id_column].astype(np.int64),
            'RA': df['RA'].astype(np.float64),
            'Dec': df['Dec'].astype(np.float64)
        })

        targets.reset_index(drop=True, inplace=True)

        pd_append_column(targets, 'prefix', prefix, 'string')

        # Add proper motion and parallax, if available
        if 'epoch' in df:
            if df['epoch'].dtype == 'string' or df['epoch'].dtype == str:
                # Convert 'J20XX.X' notation to float
                pd_append_column(targets, 'epoch', df['epoch'].apply(lambda x: float(x[1:])), np.float64)
            else:
                pd_append_column(targets, 'epoch', df['epoch'], np.float64)
        else:
            pd_append_column(targets, 'epoch', np.nan, np.float64)

        if 'pm' in df:
            pd_append_column(targets, 'pm', df['pm'], np.float64)
        else:
            pd_append_column(targets, 'pm', np.nan, np.float64)
        
        if 'pmra' in df and 'pmdec' in df:
            pd_append_column(targets, 'pmra', df['pmra'], np.float64)
            pd_append_column(targets, 'pmdec', df['pmdec'], np.float64)
        else:
            pd_append_column(targets, 'pmra', np.nan, np.float64)
            pd_append_column(targets, 'pmdec', np.nan, np.float64)

        if 'parallax' in df:
            pd_append_column(targets, 'parallax', df['parallax'], np.float64)
        else:
            pd_append_column(targets, 'parallax', np.nan, np.float64)

        if 'rv' in df:
            pd_append_column(targets, 'rv', df['rv'], np.float64)
        else:
            pd_append_column(targets, 'rv', np.nan, np.float64)
        
        # This is a per-object penalty for observing calibration targets
        if penalty is not None:
            pd_append_column(targets, 'penalty', penalty, np.int32)
        elif 'penalty' in df.columns:
            pd_append_column(targets, 'penalty', df['penalty'], np.int32)
        else:
            pd_append_column(targets, 'penalty', 0, np.int32)

        if prefix in self.__netflow_options.science_prefix:
            # Override exposure time, if specified
            if exp_time is not None:
                pd_append_column(targets, 'exp_time', exp_time, np.float64)
            else:
                pd_append_column(targets, 'exp_time', df['exp_time'], np.float64)

            # Override priority, if specified
            if priority is not None:
                pd_append_column(targets, 'priority', priority, pd.Int32Dtype())
            else:
                pd_append_column(targets, 'priority', df['priority'], pd.Int32Dtype())

            pd_append_column(targets, 'class', pd.NA, 'string')
            priority_mask = ~targets['priority'].isna()

            if prefix == 'sci':
                targets.loc[priority_mask, 'class'] = \
                    targets.loc[priority_mask, ['prefix', 'priority']].apply(
                        lambda r: f"{r['prefix']}_P{r['priority']}",
                        axis=1).astype('string')
            else:
                targets.loc[priority_mask, 'class'] = prefix
        else:
            # Calibration targets have no prescribed exposure time and priority
            pd_append_column(targets, 'exp_time', np.nan, np.float64)
            pd_append_column(targets, 'priority', pd.NA, pd.Int32Dtype())

            # Calibration targets have the same target class as the prefix
            pd_append_column(targets, 'class', prefix, 'string')

        # We allow for a non-observation cost to be specified on a target-by-target
        # basis. This is only used when the `per_target_non_obs_cost` option is turned on,
        # which is off by default.
        # If the non-observation cost is not specified on a per target basis, we
        # use the value from the configuration for the entire target class.
        if non_observation_cost is not None:
            pd_append_column(targets, 'non_observation_cost', non_observation_cost, np.int32)
        elif 'non_observation_cost' in df.columns:
            pd_append_column(targets, 'non_observation_cost', df['non_observation_cost'], np.int32)
        else:
            # The non-observation cost is not specified in the input catalog
            # In this case, we should use the non-observation cost from the configuration for each target class
            pd_append_column(targets, 'non_observation_cost', 0, np.int32)

            for target_class in targets['class'].unique():
                if target_class not in self.__netflow_options.target_classes:
                    raise NetflowException(f"Target class `{target_class}` not found in the configuration.")
                
                noc = self.__netflow_options.target_classes[target_class].non_observation_cost
                pd_update_column(targets, targets['class'] == target_class, 'non_observation_cost', noc, np.int32)
            
        # Validate the new targets
        sci_mask = targets['prefix'].isin(self.__netflow_options.science_prefix)
        if np.any(targets['priority'].isna() & sci_mask):
            raise NetflowException('Science targets must have a priority.')

        if np.any(targets['exp_time'].isna() & sci_mask):
            raise NetflowException('Science targets must have an exposure time.')

        if prefix == 'sci':
            if np.any((targets['exp_time'] <= 0) & sci_mask):
                raise NetflowException('Science targets must have a non-zero exposure time.')

        return targets
        
    def append_targets(self, targets, prefix):
        """
        Add targets from an input target list to the internal list of targets.

        Parameters
        ----------
        targets : pd.DataFrame
            A target list preprocessed by `get_catalog_data`.
        """

        # Append to the existing list
        if self.__targets is None:
            self.__targets = targets
            idx = np.array(self.__targets.index)
        else:
            n = len(self.__targets)
            self.__targets = pd.concat([self.__targets, targets], ignore_index=True)
            idx = np.array(self.__targets.index[n:])

        logger.info(f'Added {len(targets)} targets with prefix `{prefix}` to the target list.')
        
        return idx
            
    def append_science_targets(self, catalog, exp_time=None, priority=None, mask=None, filter=None, selection=None):
        """Add science targets"""

        # NOTE: This function does not cross-match the input catalog with the targets
        #       already added. Use the netflow script to achieve this.

        targets = self.get_catalog_data(catalog,
                                        prefix='sci',
                                        exp_time=exp_time,
                                        priority=priority,
                                        mask=mask,
                                        filter=filter,
                                        selection=selection)
        
        return self.append_targets(targets, prefix='sci')

    def append_sky_targets(self, sky, mask=None, filter=None, selection=None):
        """Add sky positions"""

        targets = self.get_catalog_data(sky,
                                        prefix='sky',
                                        mask=mask,
                                        filter=filter,
                                        selection=selection)
        
        return self.append_targets(targets, prefix='sky')

    def append_fluxstd_targets(self, fluxstd, mask=None, filter=None, selection=None):
        """Add flux standard positions"""

        # NOTE: This function does not cross-match the input catalog with the targets
        #       already added. Use the netflow script to achieve this.

        targets = self.get_catalog_data(fluxstd,
                                        prefix='cal',
                                        mask=mask,
                                        filter=filter,
                                        selection=selection)

        return self.append_targets(targets, prefix='cal')

    def update_targets(self, input_targets, prefix, idx2):
        """
        Update existing targets in the target list based on the contents of `input_targets`.

        This is called on targets which have already been added to the target list but appear in
        another input source list after cross-matching.

        The data frame `input_targets` must contain the same number of rows as the number
        of indexes in `idx2`.
        
        Parameters
        ==========
        input_targets: pd.DataFrame
            The catalog containing the updated target information.
        prefix: str
            The prefix to be used for the target class. This is usually `sci` for science
            targets and `cal` for calibration targets.
        idx2:
            The index of the targets in the internal target list.
        """

        targets = self.__targets

        # targetid and the coordinates are not updated

        # epoch, pm and its components can only be updated together
        # the primary columns are pm or pmra and pmdec to consider here

        update_mask = np.array(targets.loc[idx2, 'pm'].isna() &
                               (targets.loc[idx2, 'pmra'].isna() |
                                targets.loc[idx2, 'pmdec'].isna()))

        if 'epoch' in input_targets:
            if input_targets['epoch'].dtype == 'string' or input_targets['epoch'].dtype == str:
                # Convert 'J20XX.X' notation to float
                pd_update_column(targets, idx2[update_mask], 'epoch',
                                 input_targets['epoch'][update_mask].apply(
                                     lambda x: float(x[1:])),
                                 np.float64)
            else:
                pd_update_column(targets, idx2[update_mask], 'epoch',
                                 input_targets['epoch'][update_mask],
                                 np.float64)

        if 'pm' in input_targets:
            pd_update_column(targets, idx2[update_mask], 'pm', input_targets['pm'][update_mask], np.float64)
        
        if 'pmra' in input_targets and 'pmdec' in input_targets:
            pd_update_column(targets, idx2[update_mask], 'pmra', input_targets['pmra'][update_mask], np.float64)
            pd_update_column(targets, idx2[update_mask], 'pmdec', input_targets['pmdec'][update_mask], np.float64)

        if 'parallax' in input_targets:
            pd_update_column(targets, idx2[update_mask], 'parallax', input_targets['parallax'][update_mask], np.float64)

        if prefix in self.__netflow_options.science_prefix:
            # Only update exp_time if it is still NA or the new exp_time is larger than
            # the existing one
            update_mask = np.array(targets.loc[idx2, 'exp_time'].reset_index(drop=True).isna() |
                                   ((input_targets['exp_time'].reset_index(drop=True)) >
                                    (targets.loc[idx2, 'exp_time'].reset_index(drop=True))).fillna(False))
            pd_update_column(targets, idx2[update_mask], 'exp_time', input_targets['exp_time'][update_mask], np.float64)

            # Only update priority if it is still NA or the new priority is smaller than
            # the existing one
            update_mask = np.array(targets.loc[idx2, 'priority'].reset_index(drop=True).isna() |
                                   ((input_targets['priority'].reset_index(drop=True)) <
                                    (targets.loc[idx2, 'priority'].reset_index(drop=True))).fillna(False))
            pd_update_column(targets, idx2[update_mask], 'priority', input_targets['priority'][update_mask], np.int32)

            # Update the target class labels if the priority is updated.
            # When we have a priority, it must be a science target.
            # If some science targets are also flux standards, we allocate fibers to them
            # as flux standards. Care must be taken to process these as science targets by the GA pipeline.

            priority_mask = np.array(~targets.loc[idx2[update_mask], 'priority'].isna() &
                                     (targets.loc[idx2[update_mask], 'prefix'] == 'sci'))
            targets.loc[idx2[update_mask][priority_mask], 'class'] = \
                targets.loc[idx2[update_mask][priority_mask], ['prefix', 'priority']].apply(
                    # lambda r: f"{r['prefix']}_P{r['priority']}",
                    lambda r: f"sci_P{r['priority']}",
                    axis=1).astype('string')

            # TODO: implement option to re-allocate flux standards as science targets but
            #       it requires changes in how the number of targets are counted for the
            #       purposes of target class minimum and maximum, as well as cobra group minima and maxima

            # Update the non-observation cost based on the values in the catalog
            update_mask = np.array(targets.loc[idx2, 'non_observation_cost'].reset_index(drop=True).isna() |
                                    ((input_targets['non_observation_cost'].reset_index(drop=True)) >
                                    (targets.loc[idx2, 'non_observation_cost'].reset_index(drop=True))).fillna(False))
            pd_update_column(targets,
                             idx2[update_mask],
                             'non_observation_cost',
                             input_targets['non_observation_cost'][update_mask],
                             np.int32)

            # TODO: Update the number of done visits
            pass
        else:
            # Calibration targets don't override exp_time or priority
            pass

        logger.info(f'Updated {idx2.size} targets with prefix `{prefix}` in the target list.')

    def validate_targets(self):
        # TODO: Make sure cobra location/instrument group constraints are satisfiable
                
        self.__validate_science_targets()
        self.__validate_fluxstd_targets()
        self.__validate_sky_targets()

    def __validate_science_targets(self):
        # Make sure all science targets have a priority class and exposure time
        for column, ignore in zip(['priority', 'exp_time'],
                                  [self.__debug_options.ignore_missing_priority,
                                   self.__debug_options.ignore_missing_exp_time]):
            
            mask = (self.__targets['prefix'].isin(self.__netflow_options.science_prefix)) & (self.__targets[column].isna())
            if mask.any():
                msg = f"Missing {column} for {mask.sum()} science targets."
                if not ignore:
                    raise NetflowException(msg)
                else:
                    logger.warning(msg)

    def __validate_fluxstd_targets(self):
        pass

    def __validate_sky_targets(self):
        pass

    def __validate_cobra_group_limits(self):
        """
        Make sure cobra location/instrument group constraints are satisfiable
        """

        pass

    def init(self, resume=False, save=True):
        """Preprocess the targets"""
        
        logger.info("Validating netflow configuration")
        self.__validate_options()

        logger.info("Preprocessing target list")
        self.__calculate_exp_time()
        self.__calculate_target_visits()
        self.__validate_target_visits()

    def build(self, resume=False, save=True):
        """Construct the ILP problem"""

        self.__check_pointing_visibility()

        # Create a cache of target properties to optimize data access
        cache_loaded = False
        if resume:
            try:
                self.__load_target_cache()
                cache_loaded = True
            except FileNotFoundError as e:
                pass

        if not cache_loaded:
            self.__cache_targets()
            if save:
                self.__save_target_cache()

        # Forbidden targets and target pairs
        self.__forbidden_targets = self.__get_forbidden_targets_config()
        self.__forbidden_pairs = self.__get_forbidden_pairs_config()

        # Cobras positioned too close to a black dot can get a penalty
        self.__black_dots = self.__get_black_dot_config()

        # Generate the list of visits
        self.__create_visits()

        target_fp_pos_loaded = False
        if resume:
            try:
                self.__load_target_fp_pos()
                target_fp_pos_loaded = True
            except FileNotFoundError as e:
                pass

        if not target_fp_pos_loaded:
            self.__calculate_target_fp_pos()
            if save:
                self.__save_target_fp_pos()

        # Run a few sanity checks
        self.__check_target_fp_pos()

        # Calculate visibility and collisions
        visibilities_loaded = False
        if resume:
            try:
                self.__load_visibility()
                visibilities_loaded = True
            except FileNotFoundError as e:
                pass
        
        if not visibilities_loaded:
            self.__calculate_visibility()
            if save:
                self.__save_visibility()

        collisions_loaded = False
        if resume:
            try:
                self.__load_collisions()
                collisions_loaded = True
            except FileNotFoundError as e:
                pass

        if not collisions_loaded:
            self.__calculate_collisions()
            if save:
                self.__save_collisions()

        # Convert the netflow graph into a linear program
        problem_loaded = False
        if resume:
            try:
                self.load_problem()
                problem_loaded = True
            except FileNotFoundError as e:
                pass

        if not problem_loaded:
            self.__build_ilp_problem()
            if save:
                self.save_problem()

    def __check_pointing_visibility(self, airmass_limit=3.0):
        """
        Verify that all pointings are visible in the sky at obs_time and that
        the airmass is not extremely high.
        """

        for p in self.__pointings:
            # Convert pointing center into azimuth and elevation (+ instrument rotator angle)
            alt, az, inr = self.__instrument.radec_to_altaz(p.ra, p.dec, p.posang, p.obs_time)

            # Calculate airmass from AZ
            airmass = 1.0 / np.sin(np.radians(90.0 - alt))

            assert az > 0, f"Pointing is below the horizon."
            assert airmass < airmass_limit, f"Pointing is at a high airmass of {airmass:.2f}."
                    
    def __calculate_exp_time(self):
        """
        Calculate the exposure time for each object.
        """

        # Since targeting is done in fix quanta in time, make sure all pointings and visits
        # have the same exposure time

        exp_time = None
        for p in self.__pointings:
            if exp_time is None:
                exp_time = p.exp_time.to_value('s')
            elif exp_time != p.exp_time.to_value('s'):
                raise NetflowException("Exposure time of every pointing and visit should be the same.")

        if exp_time is None:
            raise NetflowException("No exposure time provided for the pointings.")

        self.__visit_exp_time = exp_time
            
    def __calculate_target_visits(self):
        """
        Calculate the number of required visits for each target.

        Assumes the same exposure time for each visit.
        """

        # TODO: The number of done visits should come from the science target lists

        sci_mask = self.__targets['prefix'].isin(self.__netflow_options.science_prefix)
        cal_mask = self.__targets['prefix'].isin(self.__netflow_options.calibration_prefix)

        # Reset the number of done visits for new targets
        if 'done_visits' not in self.__targets.columns:
            pd_append_column(self.__targets, 'done_visits', 0, np.int32)

        if 'req_visits' not in self.__targets.columns:
            pd_append_column(self.__targets, 'req_visits', 0, np.int32)

        # Some targets have already be observed
        done_mask = (self.__targets['done_visits'] > 0) & (self.__targets['done_visits'] >= self.__targets['req_visits'])

        # Calculate the number of required visits from the exposure time
        # The default is 0 for calibration targets. Some targets are both flux
        # standards and science targets, so we calculate the number of required visits
        # based on the presence of exp_time

        exp_time_mask = self.__targets['exp_time'].notna()
        mask = (~done_mask | cal_mask) & exp_time_mask
        pd_update_column(
            self.__targets, mask, 'req_visits',
            np.ceil(self.__targets['exp_time'][mask] / self.__visit_exp_time), dtype=np.int32)

        hist = np.bincount(self.__targets[sci_mask]['req_visits'])
        logger.info(f'Histogram of required visits: {hist}')

        hist = np.bincount(self.__targets[sci_mask]['done_visits'])
        logger.info(f'Histogram of done visits: {hist}')

    def __validate_target_visits(self):
        """
        Make sure that there are not targets that would require more visits then availables
        """

        # TODO: take done_visits into account

        max_req_visits = self.__targets['req_visits'].max()
        for p in self.__pointings:
            if p.nvisits < max_req_visits:
                raise NetflowException('Some science targets require more visits than provided.')

    def __cache_targets(self):
        """Extract the contents of the Pandas DataFrame for faster indexed access."""

        # Determine the target mask that filters out targets
        # that are outside the field of view
        mask = self.__filter_targets(self.__targets, self.__pointings)
        mask_any = np.any(np.stack(mask, axis=-1), axis=-1)

        self.__target_mask = mask
        self.__target_cache_mask = mask_any

        # Map target index to cache index
        self.__target_to_cache_map = np.full(len(self.__targets), -1, dtype=np.int32)
        self.__target_to_cache_map[mask_any] = np.arange(np.sum(mask_any))
        
        # Map cache index to target index
        self.__cache_to_target_map = np.where(mask_any)[0]

        self.__target_cache = SimpleNamespace(
                key = np.array(self.__targets['key'][mask_any].astype(str)),
                target_idx = np.array(self.__targets.index[mask_any].astype(np.int64)),
                id = np.array(self.__targets['targetid'][mask_any].astype(np.int64)),
                ra = np.array(self.__targets['RA'][mask_any].astype(np.float64)),
                dec = np.array(self.__targets['Dec'][mask_any].astype(np.float64)),
                pm = np.array(self.__targets['pm'][mask_any].astype(np.float64)),
                pmra = np.array(self.__targets['pmra'][mask_any].astype(np.float64)),
                pmdec = np.array(self.__targets['pmdec'][mask_any].astype(np.float64)),
                parallax = np.array(self.__targets['parallax'][mask_any].astype(np.float64)),
                epoch = np.array(self.__targets['epoch'][mask_any].astype(np.float64)),
                target_class = np.array(self.__targets['class'][mask_any].astype(str)),
                prefix = np.array(self.__targets['prefix'][mask_any].astype(str)),
                req_visits = np.array(self.__targets['req_visits'][mask_any].astype(np.int32)),
                done_visits = np.array(self.__targets['done_visits'][mask_any].astype(np.int32)),
                penalty = np.array(self.__targets['penalty'][mask_any].astype(np.int32)),
                non_observation_cost = np.array(self.__targets['non_observation_cost'][mask_any].astype(np.int32)),
            )
        
    def __get_target_cache_filename(self, filename=None):
        return self.__get_prefixed_filename(filename, 'netflow_target_cache.pkl')
        
    def __save_target_cache(self, filename=None):
        filename = self.__get_target_cache_filename(filename)
        
        with open(filename, 'wb') as f:
            data = (self.__target_cache,
                    self.__target_mask,
                    self.__target_cache_mask,
                    self.__target_to_cache_map,
                    self.__target_to_cache_map,
                    self.__cache_to_target_map)
            pickle.dump(data, f)

    def __load_target_cache(self, filename=None):
        filename = self.__get_target_cache_filename(filename)

        with open(filename, 'rb') as f:
            (self.__target_cache,
             self.__target_mask,
             self.__target_cache_mask,
             self.__target_to_cache_map,
             self.__target_to_cache_map,
             self.__cache_to_target_map) = pickle.load(f)

    def __create_visits(self):
        self.__visits = []
        visit_idx = 0
        for pointing_idx, p in enumerate(self.__pointings):
            for i in range(p.nvisits):
                # TODO: define cost to prefer early or late observation of targets
                v = Visit(visit_idx, pointing_idx, p, visit_cost=0)
                self.__visits.append(v)
                visit_idx += 1

    def __calculate_target_fp_pos(self):
        """
        Calculate the focal plane positions for each pointing. Limit the targets to
        those that are within the field of view of the PFI.
        """

        epoch = self.__get_netflow_option(self.__netflow_options.epoch, None)

        ignore_proper_motion = self.__debug_options.ignore_proper_motion
        if ignore_proper_motion:
            logger.warning('Proper motion is ignored when calculating focal plane coordinates!')

        self.__target_fp_pos = []
        self.__fp_pos_to_cache_map = []
        self.__cache_to_fp_pos_map = []

        for pidx, p in enumerate(self.__pointings):
            # Select targets that are within the field of view of the pointing
            mask = self.__target_mask[pidx][self.__target_cache_mask]

            # Create the mappings between the target cache index and the focal plane positions index
            map1 = np.where(mask)[0]
            map2 = np.full(len(mask), -1, dtype=np.int32)
            map2[map1] = np.arange(len(map1))
            self.__fp_pos_to_cache_map.append(map1)
            self.__cache_to_fp_pos_map.append(map2)
            
            # Calculate the focal plane positions of the targets for the pointing
            fp_pos = self.__calculate_fp_pos(p, 
                                             self.__target_cache.ra[mask],
                                             self.__target_cache.dec[mask],
                                             pmra=self.__target_cache.pmra[mask],
                                             pmdec=self.__target_cache.pmdec[mask],
                                             parallax=self.__target_cache.parallax[mask],
                                             epoch=epoch,
                                             ignore_proper_motion=ignore_proper_motion)
            
            self.__target_fp_pos.append(fp_pos)

    def __get_target_fp_pos_filename(self, filename=None):
        return self.__get_prefixed_filename(filename, 'netflow_target_fp_pos.pkl')

    def __save_target_fp_pos(self, filename=None):
        filename = self.__get_target_fp_pos_filename(filename)
        
        with open(filename, 'wb') as f:
            pickle.dump((self.__target_fp_pos, self.__fp_pos_to_cache_map, self.__cache_to_fp_pos_map), f)

    def __load_target_fp_pos(self, filename=None):
        filename = self.__get_target_fp_pos_filename(filename)
        
        with open(filename, 'rb') as f:
            self.__target_fp_pos, self.__fp_pos_to_cache_map, self.__cache_to_fp_pos_map = pickle.load(f)

    def __check_target_fp_pos(self):
        # Verify if each pointing contains a reasonable number of targets
        # A target is accessible if it is within the 400mm radius of the PFI

        for fp in self.__target_fp_pos:
            n = np.sum(np.abs(fp) < 400)
            assert n > 100, f"Pointing contains only {n} targets within the PFI radius."

            for m in [ np.min(fp.real), np.max(fp.real), np.min(fp.imag), np.max(fp.imag) ]:
                assert np.abs(m) > 190, f"The focal plane might not be fully covered by the target catalog."

    def __make_name_full(self, *parts):
        ### SLOW ### 4M calls!
        return "_".join(str(x) for x in parts)

    def __make_name_short(self, *parts):
        name = hex(self.__name_counter)
        self.__name_counter += 1
        return name
    
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
            
            T_i = dict(),                    # Target inflow (only science targets), key: (target_idx)
            T_o = defaultdict(list),         # Target outflows (only science targets), key: (target_idx)
            T_sink = dict(),                 # Target sinks (only science targets) key: (target_idx)

            Tv_i = dict(),                   # Target visit inflow, key: (target_idx, visit_idx)
            Tv_o = defaultdict(list),        # Target visit outflows, key: (target_idx, visit_idx)
            Cv_i = defaultdict(list),        # Cobra visit inflows, key: (cobra_idx, visit_idx)            

            CG_i = defaultdict(list),        # Cobra group inflow variables, for each visit, key: (name, ivis, gidx)
        )

    def __add_variable(self, name, lo=None, hi=None):
        if self.__use_named_variables and name in self.__variables.all:
            raise ValueError('Duplicate variable name')
        
        v = self.__problem.add_variable(name, lo=lo, hi=hi)
        self.__variables.all[name] = v
        return v
    
    def __add_variable_array(self, name, indexes, lo=None, hi=None, cost=None):
        if self.__use_named_variables and name in self.__variables.all:
            raise ValueError('Duplicate variable name')

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
            Tv_o_broken = dict(),           # Broken cobra constraints, key: (tidx, visit_idx)
            Tv_o_forb = dict(),             # Forbidden pairs, key: (tidx1, tidx2, visit_idx)

            Cv_CG_min = dict(),             # Cobra group minimum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            Cv_CG_max = dict(),             # Cobra group maximum target constraints, key: (cobra_group_name, visit_idx, cobra_group)
            Cv_i_sum = dict(),              # At most one target per cobra, key: (cobra_idx, visit_idx)
            Cv_i_max = dict(),              # Hard upper limit on the number of assigned fibers, key (visit_idx)
            Tv_i_Tv_o_sum = dict(),         # Target visit nodes must be balanced, key: (target_idx, visit_idx)
            T_i_T_o_sum = dict(),           # Inflow and outflow at every T node must be balanced, key: (target_idx)
            
            Tv_i_sum = dict(),              # Science program time budget, key: (budget_idx)
        )

    def __add_constraint(self, name, constraint):
        if self.__use_named_variables and name in self.__constraints.all:
            raise ValueError('Duplicate constraint name')
        
        self.__constraints.all[name] = constraint
        
        if isinstance(constraint, tuple):
            self.__problem.add_linear_constraint(name, *constraint)
        else:
            self.__problem.add_constraint(name, constraint)

    def __build_ilp_problem(self):
        """
        Construct the ILP problem by defining the variables and constraints.
        """

        self.__use_named_variables = self.__get_netflow_option(self.__netflow_options.use_named_variables, False)
        
        if self.__use_named_variables:
            self.__make_name = self.__make_name_full
        else:
            self.__name_counter = 0
            self.__make_name = self.__make_name_short

        # Number of visits
        nvisits = len(self.__visits)

        self.__init_problem()
        self.__init_variables()
        self.__init_constraints()

        logger.info("Creating network topology")

        # Create science target class variables which are independent
        # of visit because we want multiple visits of the same target.
        # > STC_sink_{target_class}
        self.__create_science_target_class_variables()

        # Create science target variables, for each visit
        # > STC_T[target_idx]
        # > T_sink[target_idx]
        self.__create_science_target_variables()
        
        # Create the variables for each visit of every pointing
        for visit_idx, visit in enumerate(self.__visits):
            logger.info(f'Processing visit {visit_idx + 1}.')

            # Create calibration target class variables and define cost
            # for each calibration target class. These are different for each visit
            # because we don't necessarily need the same calibration targets at each visit.
            # > CTCv_sink_{target_class}_{visit_idx}
            logger.debug("Creating calibration target class variables.")
            self.__create_calib_target_class_variables(visit)

            # > T_Tv_{visit_idx}[target_idx]
            # > CTCv_Tv_{visit_idx}[target_idx]
            # > Tv_Cv_{visit_idx}_{cobra_idx}[target_idx]
            logger.debug("Creating target and cobra visit variables.")
            self.__create_visit_variables(visit)
            self.__problem.update()

            # > Tv_o_coll_{?}_{?}_{visit_idx} (endpoint collisions)
            # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx} (elbow collisions)
            # > Tv_o_broken_{cobra_idx1}_{cobra_idx2}_{target_idx}_{visit_idx} (collision with broken cobras)
            logger.debug("Creating cobra collision constraints.")
            self.__create_cobra_collision_constraints(visit)

            logger.debug("Adding cobra non-allocation cost terms.")
            self.__add_fiber_non_allocation_cost(visit)

            # Add constraints for forbidden targets and forbidden pairs
            if not self.__debug_options.ignore_forbidden_targets:
                if self.__forbidden_targets is not None and len(self.__forbidden_targets) > 0:
                    # > Tv_o_forb_{tidx}_{tidx}_{visit_idx}
                    logger.debug("Adding forbidden target constraints")
                    self.__create_forbidden_target_constraints(visit)
            else:
                 logger.debug("Ignored forbidden target constraints")

            if not self.__debug_options.ignore_forbidden_pairs:
                if self.__forbidden_pairs is not None and len(self.__forbidden_pairs) > 0:
                    # > Tv_o_forb_{tidx1}_{tidx2}_{visit_idx}
                    logger.debug("Adding forbidden pair constraints")
                    self.__create_forbidden_pair_constraints(visit)
            else:
                logger.debug("Ignored forbidden pair constraints")
        
        logger.info("Adding constraints")

        # TODO: add log message in function beyond this point

        # The total number of science targets per target class must be balanced and
        # the minimum and maximum number of targets to be observed within
        # a science target class can be optionally constrained.
        # > STC_o_sum_{target_class}
        # > STC_o_min_{target_class}
        # > STC_o_max_{target_class}
        self.__create_science_target_class_constraints()

        for visit_idx, visit in enumerate(self.__visits):
            # Every calibration target class must be observed a minimum number of times
            # every visit
            # > CTCv_o_sum_{target_class}_{visit_idx}
            # > CTCv_o_min_{target_class}_{visit_idx}
            # > CTCv_o_max_{target_class}_{visit_idx}
            self.__create_calibration_target_class_constraints(visit_idx)

        # Inflow and outflow at every T node must be balanced
        # > T_i_T_o_sum_{target_idx}
        self.__create_science_target_constraints()

        # Inflow and outflow at every Tv node must be balanced
        # > Tv_i_Tv_o_sum_{target_idx}_{visit_idx}
        self.__create_Tv_i_Tv_o_constraints()

        # Every cobra can observe at most one target per visit 
        # > Cv_i_sum_{cobra_idx}_{visit_idx}
        self.__create_Cv_i_constraints()

        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget
        # > Tv_i_sum_{budget_name}
        self.__create_time_budget_constraints()

        # Make sure that there are enough sky targets in every Cobra location group
        # and instrument region
        # > Cv_CG_min_{cg_name}_{visit_idx}_{group_idx}
        # > Cv_CG_max_{cg_name}_{visit_idx}_{group_idx}
        self.__create_cobra_group_constraints(nvisits)

        # Make sure that enough fibers are kept unassigned, if this was requested
        # > Cv_i_max_{visit_idx}
        self.__create_unassigned_fiber_constraints(nvisits)

        # Make sure that the objective is set before the model is returned
        self.__problem.set_objective()

    def __calculate_visibility(self):

        logger.info("Calculating target visibilities")

        self.__visibility = []
        for pidx, pointing in enumerate(self.__pointings):
            logger.info(f'Calculating visibility for pointing {pidx + 1}.')

            # Get the target visibility and elbow positions for each target roughly within the pointing
            fp_pos = self.__target_fp_pos[pidx]
            fpidx_to_tidx_map = self.__fp_pos_to_cache_map[pidx]

            vis = SimpleNamespace()
            vis.targets_cobras, vis.cobras_targets = self.__calculate_visibility_and_elbow(fp_pos, fpidx_to_tidx_map)
            self.__visibility.append(vis)

            # targets_cobras keys   -> index into the target_cache
            # targets_cobras values -> (cidx, elbow position)
            
            # cobras_targets keys   -> cobra index
            # cobras_targets values -> (tidx, elbow position), where tidx is the index into the target_cache

    def __get_visibility_filename(self, filename=None):
        return self.__get_prefixed_filename(filename, 'netflow_visibilities.pkl')

    def __save_visibility(self, filename=None):
        filename = self.__get_visibility_filename(filename)
        
        # Save the visibility data using pickle
        with open(filename, 'wb') as f:
            pickle.dump(self.__visibility, f)

    def __load_visibility(self, filename=None):
        filename = self.__get_visibility_filename(filename)
        
        # Load the visibility data using pickle
        with open(filename, 'rb') as f:
            self.__visibility = pickle.load(f)

    def __calculate_collisions(self):
        logger.info("Calculating cobra collisions")

        collision_distance = self.__get_netflow_option(self.__netflow_options.collision_distance, 0.0)

        self.__collisions = []
        for pidx, pointing in enumerate(self.__pointings):
            logger.info(f'Calculating collisions for pointing {pidx + 1}.')

            # Collect endpoint and elbow collisions for each target pair
            col = SimpleNamespace()
            col.endpoints = self.__get_colliding_endpoints(pidx, collision_distance)
            col.elbows = self.__get_colliding_elbows(pidx, collision_distance)
            
            # TODO
            # col.trajectories = self.__get_colliding_trajectories(pidx, collision_distance)

            # Collect endpoint and elbow collisions for targets and broken cobras
            col.broken_cobras = self.__get_colliding_broken_cobras(pidx, collision_distance)
    
            self.__collisions.append(col)

            logger.info(f'Found {len(col.endpoints)} colliding fibers and {len(col.elbows)} elbows for pointing {pidx + 1}.')
            logger.info(f'Found {len(col.broken_cobras)} colliding fibers with broken cobras for pointing {pidx + 1}.')

    def __get_collisions_filename(self, filename=None):
        return self.__get_prefixed_filename(filename, 'netflow_collisions.pkl')

    def __save_collisions(self, filename=None):
        filename = self.__get_collisions_filename(filename)

        with open(filename, 'wb') as f:
            pickle.dump(self.__collisions, f)

    def __load_collisions(self, filename=None):        
        filename = self.__get_collisions_filename(filename)

        with open(filename, 'rb') as f:
            self.__collisions = pickle.load(f)

    #region Variables and costs
        
    def __create_science_target_class_variables(self):
        for target_class in self.__netflow_options.target_classes.keys():
            if self.__netflow_options.target_classes[target_class].prefix in self.__netflow_options.science_prefix:
                # Science Target class node to sink
                self.__create_STC_sink(target_class)

    def __create_STC_sink(self, target_class):
        # Non-observed science target add to the cost with a coefficient defined
        # for the particular science target class. Higher priority targets should add
        # more to the cost to prefer targeting them.
        f = self.__add_variable(self.__make_name("STC_sink", target_class), 0, None)
        self.__variables.STC_sink[target_class] = f

        # If the non-observation cost is handled on the target class level,
        # define it now on the sink edge.
        if not self.__netflow_options.per_target_non_obs_cost:
            self.__add_cost(f * self.__netflow_options.target_classes[target_class].non_observation_cost)
        
    def __create_science_target_variables(self):
        """
        Create the variables corresponding to the STC -> T and STC -> STC_sink edges.
        These edges are created regardless of the visibility of the targets during the visits.
        """

        # Select only science targets and create the variables in batch mode
        # Ignore any science target with unknown priority class (how can this happen? this means a wrong config)
        mask = (np.isin(self.__target_cache.prefix, self.__netflow_options.science_prefix) &
                np.isin(self.__target_cache.target_class, list(self.__netflow_options.target_classes.keys())))
        
        tidx = np.where(mask)[0]
        target_class = self.__target_cache.target_class[mask]
        target_penalty = self.__target_cache.penalty[mask]
        target_non_observation_cost = self.__target_cache.non_observation_cost[mask]

        # Science Target class node to target node: STC -> T
        self.__create_STC_T(tidx, target_class, target_penalty, target_non_observation_cost)

        # Science Target node to sink: T -> T_sink
        # > T_sink[target_idx]
        self.__create_T_sink(tidx, target_class)

    def __create_STC_T(self, tidx, target_class, target_penalty, target_non_observation_cost):
        vars =  self.__add_variable_array('STC_T', tidx, 0, 1)
        for i in range(len(tidx)):
            f = vars[tidx[i]]
            self.__variables.T_i[tidx[i]] = f
            self.__variables.STC_o[target_class[i]].append(f)
            
            # TODO: figure out how to add these cost terms in batch mode
            # use something like vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

            # The observation of a target can be penalized
            if target_penalty[i] != 0:
                self.__add_cost(f * target_penalty[i])

            # If the non-observation cost is handled on a per target basis,
            # add the cost term here.
            if self.__netflow_options.per_target_non_obs_cost:
                if target_non_observation_cost[i] != 0:
                    self.__add_cost((1 - f) * target_non_observation_cost[i])
        
    def __create_T_sink(self, tidx, target_class):
        """
        Create the sink variable which is responsible for draining the flow from the target nodes
        which get less visits than required
        """

        max_visits = len(self.__visits)
        T_sink =  self.__add_variable_array('T_sink', tidx, 0, max_visits)

        for i in range(len(tidx)):
            self.__variables.T_sink[tidx[i]] = T_sink[tidx[i]]

            # TODO: this could be added at once as a LinExpr
            self.__add_cost(T_sink[tidx[i]] * self.__netflow_options.target_classes[target_class[i]].partial_observation_cost)

    def __create_calib_target_class_variables(self, visit):
        # Calibration target class visit outflows, key: (target_class, visit_idx)
        for target_class, options in self.__netflow_options.target_classes.items():
            if options.prefix in self.__netflow_options.calibration_prefix:
                f = self.__add_variable(self.__make_name("CTCv_sink", target_class, visit.visit_idx), 0, None)
                self.__variables.CTCv_sink[(target_class, visit.visit_idx)] = f

                # The non-observation cost is always handled on the
                # target class level for calibration targets because there is no
                # single variable that can be used to define the cost.
                self.__add_cost(f * options.non_observation_cost)
    
    def __create_visit_variables(self, visit):
        """
        Create the variables for the visible variables for a visit.
        """

        visibility = self.__visibility[visit.pointing_idx]
        tidx = np.array(list(visibility.targets_cobras.keys()), dtype=int)

        # Create science targets T_Tv edges in batch
        # > T_Tv_{visit_idx}[target_idx]
        sci_mask = (np.isin(self.__target_cache.prefix[tidx], self.__netflow_options.science_prefix) &
                    np.isin(self.__target_cache.target_class[tidx], list(self.__netflow_options.target_classes.keys())))
        self.__create_T_Tv(visit, tidx[sci_mask])

        # Create calibration targets CTCv_Tv edges
        # > CTCv_Tv_{visit_idx}[target_idx]
        cal_mask = np.isin(self.__target_cache.prefix[tidx], self.__netflow_options.calibration_prefix)
        self.__create_CTCv_Tv(visit, tidx[cal_mask])

        # For each Cv, generate the incoming Tv_Cv edges
        # These differ in number but can be vectorized for each cobra
        # > Tv_Cv_{visit_idx}_{cobra_idx}[target_idx]
        for cidx, tidx_elbow in visibility.cobras_targets.items():
            tidx = np.array([ ti for ti, _, _ in tidx_elbow ], dtype=int)
            self.__create_Tv_Cv_CG(visit, cidx, tidx)

    def __create_T_Tv(self, visit, tidx):
        """
        Create the T -> Tv variables for a visit. This function is called for the visible
        targets only.

        Parameters:
        -----------
        visit : Visit
            The visit object
        tidx : array
            The target indices visible by the cobras in the visit
        """

        vars = self.__add_variable_array(self.__make_name('T_Tv', visit.visit_idx), tidx, 0, 1)
        for ti in tidx:
            f = vars[ti]
            self.__variables.T_o[ti].append(f)
            self.__variables.Tv_i[(ti, visit.visit_idx)] = f
            
    def __create_CTCv_Tv(self, visit, tidx):
        cost = self.__target_cache.penalty[tidx]
        name = self.__make_name('CTCv_Tv', visit.visit_idx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]
            self.__variables.Tv_i[(ti, visit.visit_idx)] = f
            self.__variables.CTCv_o[(target_class, visit.visit_idx)].append(f)

    def __create_Tv_Cv_CG(self, visit, cidx, tidx):
        """
        Create the Tv -> Cv variables for a given visit and a cobra and a list of accessible targets.

        Parameters:
        -----------
        visit : Visit
            The visit object
        cidx : int
            The cobra index
        tidx : array
            The target indices visible by the cobra
        """
        
        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[visit.pointing_idx]
        
        # Calculate the cost for each target - cobra assignment
        cost = np.zeros_like(tidx, dtype=float)
        
        # Cost of a single visit
        if visit.visit_cost is not None:
            cost += visit.visit_cost

        # Focal plane positions for the current visit
        fp_pos = self.__target_fp_pos[visit.pointing_idx][tidx_to_fpidx_map[tidx]]

        # Cost of moving the cobra away from the center
        cobra_move_cost = self.__netflow_options.cobra_move_cost
        if cobra_move_cost is not None:
            dist = self.__get_cobra_center_distance(cidx, fp_pos)
            cost += cobra_move_cost(dist)
        
        # Cost of closest black dots for each cobra
        if self.__black_dots.black_dot_penalty is not None:
            dist = self.__get_closest_black_dot_distance(cidx, fp_pos)
            cost += self.__black_dots.black_dot_penalty(dist)

        # Create LP variables
        name = self.__make_name('Tv_Cv', visit.visit_idx, cidx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        # For each accessible target, save the Tv -> Cv variable to the respective lists
        # TODO: can these be created in batch mode?
        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]

            # TODO: consider creating a separate list for each visit
            self.__variables.Cv_i[(cidx, visit.visit_idx)].append((f, ti))
            self.__variables.Tv_o[(ti, visit.visit_idx)].append((f, cidx))

            # Save the variable to the list of each cobra group to which it's relevant
            for cg_name, options in self.__netflow_options.cobra_groups.items():
                if target_class in options.target_classes:
                    self.__variables.CG_i[(cg_name, visit.visit_idx, options.groups[cidx])].append(f)

    #endregion
    #region Special cost term
            
    def __add_fiber_non_allocation_cost(self, visit):
        # If requested, penalize non-allocated fibers
        # Sum up the all Cv_i edges for the current visit and penalize its difference from
        # the total number of cobras
        fiber_non_allocation_cost = self.__get_netflow_option(self.__netflow_options.fiber_non_allocation_cost, 0)
        if fiber_non_allocation_cost != 0:
            # TODO: consider storing Cv_i organized by visit instead of cobra as well
            relevant_vars = [ var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == visit.visit_idx ]
            relevant_vars = [ var for sublist in relevant_vars for (var, ti) in sublist ]
            self.__add_cost(fiber_non_allocation_cost *
                            (self.__bench.cobras.nCobras - self.__problem.sum(relevant_vars)))
            
    #endregion
    #region Constraints

    def __create_cobra_collision_constraints(self, visit):
        ignore_endpoint_collisions = self.__debug_options.ignore_endpoint_collisions
        ignore_elbow_collisions = self.__debug_options.ignore_elbow_collisions
        ignore_broken_cobra_collisions = self.__debug_options.ignore_broken_cobra_collisions
        
        # Add constraints 
        logger.debug(f"Adding constraints for visit {visit.visit_idx + 1}")

        # Avoid endpoint or elbow collisions
        collision_distance = self.__get_netflow_option(self.__netflow_options.collision_distance, 0.0)
        
        if collision_distance > 0.0:
            if not ignore_endpoint_collisions:
                # > Tv_o_coll_{target_idx1}_{target_idx2}_{visit_idx}
                logger.debug("Adding endpoint collision constraints")
                self.__create_endpoint_collision_constraints(visit)
            else:
                logger.debug("Ignoring endpoint collision constraints")

            if not ignore_elbow_collisions:
                # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx}
                logger.debug("Adding elbow collision constraints")
                self.__create_elbow_collision_constraints(visit)
            else:
                logger.debug("Ignoring elbow collision constraints")

            if not ignore_broken_cobra_collisions:
                # > Tv_o_broken_{cobra_idx1}_{cobra_idx2}_{target_idx}_{visit_idx}
                logger.debug('Adding broken cobra collision constraints')
                self.__create_broken_cobra_collision_constraints(visit)
            else:
                logger.debug('Ignoring broken cobra collision constraints')
                
    def __create_endpoint_collision_constraints(self, visit):
        """
        Create constraints to prevent allocating fibers to targets which would cause collisions with the
        endpoints of other cobras.

        The constraint prevents allocating fibers to both targets at the same time and not by
        completely excluding the targets. Separate constraints are created for each visit to allow for more
        scheduling flexibility.
        """

        vidx = visit.visit_idx
        collisions = self.__collisions[visit.pointing_idx]
        
        for tidx1, tidx2 in collisions.endpoints:
            vars = [ v for v, cidx in self.__variables.Tv_o[(tidx1, vidx)] ] + \
                   [ v for v, cidx in self.__variables.Tv_o[(tidx2, vidx)] ]
            
            name = self.__make_name("Tv_o_coll", tidx1, tidx2, vidx)
            # constr = self.__problem.sum(vars) <= 1
            constr = ([1] * len(vars), vars, '<=', 1)
            self.__constraints.Tv_o_coll[(tidx1, tidx2, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_elbow_collision_constraints(self, visit):
        """
        Create constraints to prevent collisions between the tip of one cobra with the phi arm
        of another.

        The constraint prevents allocating fibers to both targets at the same time and not by
        completely excluding the targets. Separate constraints are created for each visit to allow for more
        scheduling flexibility.
        """

        vidx = visit.visit_idx
        collisions = self.__collisions[visit.pointing_idx]

        for (cidx1, tidx1), cidx2_tidx2_list in collisions.elbows.items():

            # Find the Tv_o variable that corresponds to tidx1 and cidx1 for visit vidx
            vars = [ f for f, ci in self.__variables.Tv_o[(tidx1, vidx)] if ci == cidx1]

            # Loop over targets seen by neighboring cobras that would cause an
            # elbow collision with cidx1
            for cidx2, tidx2 in cidx2_tidx2_list:
                
                # Find the Tv_o variables of all targets seen by neighboring cobras that would
                # cause an elbow collision with cidx1 poiting at tidx1
                vars += [ f for f, ci2 in self.__variables.Tv_o[(tidx2, vidx)] if ci2 != cidx1 ]
    
                name = self.__make_name("Tv_o_coll", tidx1, cidx1, tidx2, cidx2, vidx)
                # constr = self.__problem.sum(vars) <= 1
                constr = ([1] * len(vars), vars, '<=', 1)
                self.__constraints.Tv_o_coll[(tidx1, cidx1, tidx2, cidx2, vidx)] = constr
                self.__add_constraint(name, constr)

    def __create_broken_cobra_collision_constraints(self, visit):
        """
        Create constraints that prevent allocating fibers to targets which would cause collisions with the
        endpoints of broken cobras.

        The constraint prevents allocating fibers to these targets, because the broken cobras cannot
        be moved out of the way.
        """

        vidx = visit.visit_idx
        collisions = self.__collisions[visit.pointing_idx]
        
        # cidx is the broken cobra, ncidx is its neighbor
        # tidx_list is the list of targets visible by both that would cause a collision
        # due to cidx parked at a dead position
        # None of these can be observed which is expressed by the constraint that the
        # sum of the Tv_Cv variables must be zero.
        for (cidx, ncidx), tidx_list in collisions.broken_cobras.items():
            for tidx in tidx_list:
                vars = [ v for v, ci in self.__variables.Tv_o[(tidx, vidx)] ]
                
                name = self.__make_name("Tv_o_broken", cidx, ncidx, tidx, vidx)
                # constr = self.__problem.sum(vars) = 0
                constr = ([1] * len(vars), vars, '=', 0)
                self.__constraints.Tv_o_broken[(cidx, ncidx, tidx, vidx)] = constr
                self.__add_constraint(name, constr)

    def __create_forbidden_target_constraints(self, visit):
        """
        Create the constraints prohibiting individual targets being observed in any visit.
        """

        # All edges of visible targets, relevant for this visit
        tidx_set = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == visit.visit_idx)
        for tidx in self.__forbidden_targets:
            if tidx in tidx_set:
                flows = [ v for v, cidx in self.__variables.Tv_o[(tidx, visit.visit_idx)] ]
                name = self.__make_name("Tv_o_forb", tidx, tidx, visit.visit_idx)
                # constr = self.__problem.sum(flows) == 0
                constr = ([1] * len(flows), flows, '==', 0)
                self.__constraints.Tv_o_forb[(tidx, tidx, visit.visit_idx)] = constr
                self.__add_constraint(name, constr)

    def __create_forbidden_pair_constraints(self, visit):
        """
        Create the constraints prohibiting two targets (or individual targets) being observed
        in the same visit.
         
        This constraint is independent of the cobra collision constraints.
        """

        # All edges of visible targets, relevant for this visit
        tidx_set = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == visit.visit_idx)
        for p in self.__forbidden_pairs:
            [tidx1, tidx2] = p
            if tidx1 in tidx_set and tidx2 in tidx_set:
                flows = [ v for v, cidx in self.__variables.Tv_o[(tidx1, visit.visit_idx)] ] + \
                        [ v for v, cidx in self.__variables.Tv_o[(tidx2, visit.visit_idx)] ]
                name = self.__make_name("Tv_o_forb", tidx1, tidx2, visit.visit_idx)
                # constr = self.__problem.sum(flows) <= 1
                constr = ([1] * len(flows), flows, '<=', 1)
                self.__constraints.Tv_o_forb[(tidx1, tidx2, visit.visit_idx)] = constr
                self.__add_constraint(name, constr)

    def __create_science_target_class_constraints(self):
        
        # TODO: other than balancing the total number of targets in the class,
        #       the class minimum and maximum constraints don't make too much sense
        #       when netflow is run in stages because the limits will apply to
        #       a single stage only.
        
        # Science targets must be either observed or go to the sink
        # If a maximum on the science target is set fo the target class, enforce that
        # in a separate constraint
        # It defaults to the total number of targets but can be overridden for each target class
        ignore_science_target_class_minimum = self.__debug_options.ignore_science_target_class_minimum
        ignore_science_target_class_maximum = self.__debug_options.ignore_science_target_class_maximum
        
        for target_class in self.__netflow_options.target_classes:
            if self.__netflow_options.target_classes[target_class].prefix in self.__netflow_options.science_prefix:
                vars = self.__variables.STC_o[target_class]
                sink = self.__variables.STC_sink[target_class]
                num_targets = len(vars)

                # All STC_o edges must be balanced. We don't handle the incoming edges in the model
                # so simply the sum of outgoing edges must add up to the number of targets.
                name = self.__make_name("STC_o_sum", target_class)
                # constr = self.__problem.sum(vars + [ sink ]) == num_targets
                constr = ([1] * (len(vars) + 1), vars + [ sink ], '==', num_targets)
                self.__constraints.STC_o_sum[target_class] = constr
                self.__add_constraint(name, constr)

                # If a minimum or maximum constraint is set on the number observed targets in this
                # class, enforce it through a maximum constraint on the sum of the outgoing edges not
                # including the sink
                min_targets = self.__netflow_options.target_classes[target_class].min_targets
                if not ignore_science_target_class_minimum and min_targets is not None:
                    name = self.__make_name("STC_o_min", target_class)
                    # constr = self.__problem.sum(vars) >= min_targets
                    constr = ([1] * len(vars), vars, '>=', min_targets)
                    self.__constraints.STC_o_min[target_class] = constr
                    self.__add_constraint(name, constr)

                max_targets = self.__netflow_options.target_classes[target_class].max_targets
                if not ignore_science_target_class_maximum and max_targets is not None:
                    name = self.__make_name("STC_o_max", target_class)
                    # constr = self.__problem.sum(vars) <= max_targets
                    constr = ([1] * len(vars), vars, '<=', max_targets)
                    self.__constraints.STC_o_max[target_class] = constr
                    self.__add_constraint(name, constr)

    def __create_calibration_target_class_constraints(self, vidx):
        """
        Create constrains for calibration target classes that set the minimum and
        maximum number of required targets for each calibration target class for
        each visit.
        """
        
        ignore_calib_target_class_minimum = self.__debug_options.ignore_calib_target_class_minimum
        ignore_calib_target_class_maximum = self.__debug_options.ignore_calib_target_class_maximum
        
        for target_class in self.__netflow_options.target_classes:
            if self.__netflow_options.target_classes[target_class].prefix in self.__netflow_options.calibration_prefix:
                vars = self.__variables.CTCv_o[(target_class, vidx)]
                sink = self.__variables.CTCv_sink[(target_class, vidx)]
                num_targets = len(vars)

                # The sum of all outgoing edges must be equal the number of calibration targets within each class
                name = self.__make_name("CTCv_o_sum", target_class, vidx)
                # constr = self.__problem.sum(vars + [ sink ]) == num_targets
                constr = ([1] * (len(vars) + 1), vars + [ sink ], '==', num_targets)
                self.__constraints.CTCv_o_sum[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)

                # Every calibration target class must be observed a minimum number of times every visit
                min_targets = self.__netflow_options.target_classes[target_class].min_targets
                if not ignore_calib_target_class_minimum and min_targets is not None:
                    name = self.__make_name("CTCv_o_min", target_class, vidx)
                    # constr = self.__problem.sum(vars) >= min_targets
                    constr = ([1] * len(vars), vars, '>=', min_targets)
                    self.__constraints.CTCv_o_min[(target_class, vidx)] = constr
                    self.__add_constraint(name, constr)

                # Any calibration target class cannot be observed more tha a maximum number of times every visit
                max_targets = self.__netflow_options.target_classes[target_class].max_targets
                if not ignore_calib_target_class_maximum and max_targets is not None:
                    name = self.__make_name("CTCv_o_max", target_class, vidx)
                    # constr = self.__problem.sum(vars) <= max_targets
                    constr = ([1] * len(vars), vars, '<=', max_targets)
                    self.__constraints.CTCv_o_max[(target_class, vidx)] = constr
                    self.__add_constraint(name, constr)

    def __create_science_target_constraints(self):
        """
        Create the constraints that make sure the targets get as many visits as required and if not,
        they get penalized. The penalty is added in `__create_T_sink`.

        Optionally, allow for more visits than required. In this case we create two constraints:
        - One that requires at least the required number of visits
        - One that allows for more visits than required but sets a cap at the maximum number of visits.
        Note, that if `T_i` in zero, the target is not observed at all.
        """

        # TODO: handle already observed targets here

        allow_more_visits = self.__get_netflow_option(self.__netflow_options.allow_more_visits, False)

        for tidx, T_o in self.__variables.T_o.items():
            T_i = self.__variables.T_i[tidx]
            T_sink = self.__variables.T_sink[tidx]
            req_visits = self.__target_cache.req_visits[tidx]
            done_visits = self.__target_cache.done_visits[tidx]
            max_visits = len(self.__visits)

            if not allow_more_visits:
                # Require an exact number of visits
                name = self.__make_name("T_i_T_o_sum", tidx)
                # constr = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) == done_visits
                constr = ([ req_visits ] + [ -1 for _ in T_o ] + [ -1 ],
                          [ T_i ] + T_o + [ T_sink ], '==', done_visits)
                self.__constraints.T_i_T_o_sum[tidx] = constr
                self.__add_constraint(name, constr)
            else:
                # Allow for a larger number of visits than required
                name0 = self.__make_name("T_i_T_o_sum_0", tidx)
                name1 = self.__make_name("T_i_T_o_sum_1", tidx)

                # The number of visits (together with the number of already done visits) must be
                # at least the number of required visits
                # constr0 = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) <= done_visits
                constr0 = ([ req_visits ] + [ -1 for _ in T_o ] + [ -1 ],
                           [ T_i ] + T_o + [ T_sink ], '<=', done_visits)

                # The total number of outgoing edges must not be larger, together with the sink,
                # than the number of visits. This is true regardless how many visits have been done so far,
                # so done_visits doesn't play a role here.
                # constr1 = self.__problem.sum([ max_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) >= 0
                constr1 = ([ max_visits ] + [ -1 for _ in T_o ] + [ -1 ],
                           [ T_i ] + T_o + [ T_sink ], '>=', 0)

                self.__constraints.T_i_T_o_sum[(tidx, 0)] = constr0
                self.__constraints.T_i_T_o_sum[(tidx, 1)] = constr1
                self.__add_constraint(name0, constr0)
                self.__add_constraint(name1, constr1)

    def __create_Tv_i_Tv_o_constraints(self):
        # Inflow and outflow at every Tv node must be balanced
        for (tidx, vidx), in_var in self.__variables.Tv_i.items():
            out_vars = self.__variables.Tv_o[(tidx, vidx)]

            name = self.__make_name("Tv_i_Tv_o_sum", tidx, vidx)
            # constr = self.__problem.sum([ v for v in in_vars ] + [ -v[0] for v in out_vars ]) == 0
            constr = ([ 1 ] + [ -1 ] * len(out_vars), [ in_var ] + [ v for v, _ in out_vars ], '==', 0)
            self.__constraints.Tv_i_Tv_o_sum[(tidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_Cv_i_constraints(self):
        # Every cobra can observe at most one target per visit
        for (cidx, vidx), vars in self.__variables.Cv_i.items():
            vars = [ var for (var, tidx) in vars ]
            name = self.__make_name("Cv_i_sum", cidx, vidx)
            # constr = self.__problem.sum([ f for f in vars ]) <= 1
            constr = ([1] * len(vars), vars, '<=', 1)
            self.__constraints.Cv_i_sum[(cidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_time_budget_constraints(self):
        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget

        if not self.__debug_options.ignore_time_budget and self.__netflow_options.time_budgets is not None:
            for budget_name, options in self.__netflow_options.time_budgets.items():
                budget_target_classes = set(options.target_classes)
                budget_variables = []
                # Collect all visits of targets that fall into to budget
                for (tidx, vidx), v in self.__variables.Tv_i.items():
                    target_class = self.__target_cache.target_class[tidx]
                    if target_class in budget_target_classes:
                        budget_variables.append(v)

                # TODO: This assumes a single, per visit exposure time which is fine for now
                #       but the logic could be extended further
                # TODO: add done_visits
                name = self.__make_name("Tv_i_sum", budget_name)
                # constr = self.__visit_exp_time * self.__problem.sum([ v for v in budget_variables ]) <= 3600 * options.budget
                constr = ([self.__visit_exp_time.value] * len(budget_variables), budget_variables, '<=', 3600 * options.budget)
                self.__constraints.Tv_i_sum[budget_name] = constr
                self.__add_constraint(name, constr)

    def __create_cobra_group_constraints(self, nvisits):
        # Make sure that there are enough targets in every cobra group for each visit

        ignore_cobra_group_minimum = self.__debug_options.ignore_cobra_group_minimum
        ignore_cobra_group_maximum = self.__debug_options.ignore_cobra_group_maximum

        for cg_name, options in self.__netflow_options.cobra_groups.items():
            need_min = not ignore_cobra_group_minimum and options.min_targets is not None
            need_max = not ignore_cobra_group_maximum and options.max_targets is not None
            if need_min or need_max:                
                for vidx in range(nvisits):
                    for gidx in np.unique(options.groups):
                        variables = self.__variables.CG_i[(cg_name, vidx, gidx)]

                        if need_min and len(variables) < options.min_targets:
                            raise NetflowException(f"Insufficient number of targets in `{cg_name}` for visit {vidx} in cobra group {gidx}."
                                                   f" Expected at least {options.min_targets}, found {len(variables)}.")
                        elif need_min and options.min_targets > 0:
                            if len(variables) < 3 * options.min_targets:
                                logger.warning(f"Number of targets is very low for `{cg_name}` for visit {vidx} in cobra group {gidx}."
                                               f" Expected at least {options.min_targets}, found {len(variables)}.")

                            name = self.__make_name("Cv_CG_min", cg_name, vidx, gidx)
                            # constr = self.__problem.sum([ v for v in variables ]) >= options.min_targets
                            constr = ([1] * len(variables), variables, '>=', options.min_targets)
                            self.__constraints.Cv_CG_min[cg_name, vidx, gidx] = constr
                            self.__add_constraint(name, constr)                          

                        if need_max:
                            name = self.__make_name("Cv_CG_max", cg_name, vidx, gidx)
                            # constr = self.__problem.sum([ v for v in variables ]) <= options.max_targets
                            constr = ([1] * len(variables), variables, '<=', options.max_targets)
                            self.__constraints.Cv_CG_max[cg_name, vidx, gidx] = constr
                            self.__add_constraint(name, constr)


    def __create_unassigned_fiber_constraints(self, nvisits):
        # Make sure that enough fibers are kept unassigned, if this was requested
        # This is done by setting an upper limit on the sum of Cv_i edges

        ignore_reserved_fibers = self.__debug_options.ignore_reserved_fibers
        num_reserved_fibers = self.__get_netflow_option(self.__netflow_options.num_reserved_fibers, 0)

        if not ignore_reserved_fibers and num_reserved_fibers > 0:
            max_assigned_fibers = self.__bench.cobras.nCobras - num_reserved_fibers
            for vidx in range(nvisits):
                variables = [ vars for ((ci, vi), vars) in self.__variables.Cv_i.items() if vi == vidx ]
                variables = [ var for sublist in variables for (var, tidx) in sublist ]

                name = self.__make_name("Cv_i_max", vidx)
                # constr = self.__problem.sum(variables) <= max_assigned_fibers
                constr = ([1] * len(variables), variables, '<=', max_assigned_fibers)
                self.__constraints.Cv_i_max[vidx] = constr
                self.__add_constraint(name, constr)

    #endregion

    def solve(self):
        if self.__problem.solve():
            # If the processing is resumed and the variables are empty we
            # need to restore them
            if self.__variables is None or self.__constraints is None:
                self.__restore_variables()
                self.__restore_constraints()
        else:
            if len(self.__problem.infeasible_constraints) > 0:
                # TODO: Find the infeasible constraints and give detailed information about them
                pass

            raise NetflowException("Failed to solve the netflow problem.")

    def extract_assignments(self):
        """
        Extract the fiber assignments from an LP solution by converting variable values
        into meaningful fiber assignment info.
        
        The assignments are stored in two lists:
        - target_assignments: a list of dictionaries for each visit that maps target indices
          to cobra indices
        - cobra_assignments: a list of arrays for each visit that contain an array with the id
          of associated the targets
        """

        nvisits = len(self.__visits)
        ncobras = self.__bench.cobras.nCobras

        self.__target_assignments = [ {} for _ in range(nvisits) ]
        self.__cobra_assignments = [ np.full(ncobras, -1, dtype=int) for _ in range(nvisits) ]

        def set_assigment(tidx, cidx, vidx):
            # Make sure we don't have any duplicate assignments
            assert tidx not in self.__target_assignments[vidx]
            assert self.__cobra_assignments[vidx][cidx] == -1

            self.__target_assignments[vidx][tidx] = cidx
            self.__cobra_assignments[vidx][cidx] = tidx

        if self.__variables is not None:
            # This only works when the netflow problem is fully built
            for (tidx, vidx), vars in self.__variables.Tv_o.items():
                for (f, cidx) in vars:
                    if self.__problem.is_one(f):
                        set_assigment(tidx, cidx, vidx)
        else:
            # This also works when only the LP problem is loaded back from a file
            # but it's much slower. It also requires that the problem is build with
            # full variable names.
            for k1, v1 in self.__problem._variables.items():
              if k1.startswith("Tv_Cv_"):
                    if self.__problem.is_one(v1) > 0:
                        _, _, tidx, cidx, vidx = k1.split("_")
                        set_assigment(int(tidx), int(cidx), int(vidx))

        # For each visit, check if there are any duplicate assignments
        # __target_assignments contains a dict of (tidx, cidx) pairs but
        # both must be unique for a given vidx
        for vidx in range(len(self.__target_assignments)):
            # Create a list of fiber indices, then sort them and make sure there's no duplicates
            assert not np.any(np.diff(np.sort([ cidx for tidx, cidx in self.__target_assignments[vidx].items() ])) == 0)

    def update_done_visits(self):
        """
        Update the number of done visits in the targets data frame.
        """

        # self.__target_assignments is keyed by the target_idx and contains the cobra_idx (both 0-based)

        for vidx, visit in enumerate(self.__visits):
            # Map the indices of target assignments to the target list index
            tidx = np.array([ ti for ti in self.__target_assignments[vidx].keys() ])
            tidx = self.__cache_to_target_map[tidx]

            # Update the number of done visits
            self.__targets.loc[tidx, 'done_visits'] = self.__targets.loc[tidx, 'done_visits'] + 1

    def simulate_collisions(self):
        """
        Run a cobra collision simulation for the current fiber assignments.
        """

        if self.__cobra_assignments is None:
            raise NetflowException("No fiber assignments have been calculated yet. Run Netflow.solve() first.")
        
        logger.info("Simulating trajectory collisions...")

        self.__cobra_collisions = {}

        for vidx, visit in enumerate(self.__visits):
            # Collect fiber positions ordered by cobra id
            fp_pos = np.full(self.__bench.cobras.nCobras, NULL_TARGET_POSITION)

            # Assume the position of broken cobras from the calibration model
            # These include cobras that are broken or have a broken fibers
            # mask = self.__instrument.bench.cobras.hasProblem
            # fp_pos[mask] = self.__instrument.bench.cobras.home0[mask]

            # Assigned cobras
            mask = self.__cobra_assignments[vidx] != -1
            fp_pos[mask] = self.__target_fp_pos[visit.pointing_idx][self.__cache_to_fp_pos_map[visit.pointing_idx][self.__cobra_assignments[vidx][mask]]]

            targets = TargetGroup(fp_pos, ids=None, priorities=None)

            # simulator = CollisionSimulator(self.__instrument, targets)
            # simulator.run()
            # trajectory_collisions = simulator.trajectory_collisions
            
            
            simulator = CollisionSimulator2(self.__instrument.bench, self.__instrument.cobra_coach, targets)
            simulator.run()
            trajectory_collisions = self.__instrument.bench.cobraAssociations[:, simulator.associationCollisions]
            endpoint_collisions = self.__instrument.bench.cobraAssociations[:, simulator.associationEndPointCollisions]

            logger.warning(f'{trajectory_collisions.shape[1]} trajectory collisions found for visit {vidx}.')
            logger.info(f'The colliding cobras are {trajectory_collisions}')
            
            # Extract the collisions
            self.__cobra_collisions[vidx] = trajectory_collisions

            # Verify if the collisions are correct

            # Some of the target ids are -1 which means the cobra is unassigned but
            # still collides with another cobra. This is likely because the cobra is
            # broken. Verify if this is the case.
            cobra_assignments = self.__cobra_assignments[vidx]
            tidx = cobra_assignments[trajectory_collisions]
            cidx = trajectory_collisions[tidx == -1]
            good_cobra = ~self.__instrument.bench.cobras.hasProblem[cidx]
            if not np.any(good_cobra):
                logger.warning(f'Collision detected for working cobras that are not assigned to a target.'
                               f'Visit index {vidx}, cobra ID {cidx[good_cobra] + 1}.')

    def unassign_colliding_cobras(self):
        """
        When a cobra collision is detected, unassign the cobra with the lower priority target. It
        will cause sending it to home position.
        """

        # TODO: we could do better collision resolution strategies but this will do for now

        logger.info("Unassigning cobras causing collisions.")

        if self.__cobra_collisions is None:
            raise NetflowException("No cobra collisions have been calculated yet. Run Netflow.simulate_collisions() first.")
        
        for vidx, visit in enumerate(self.__visits):
            collisions = self.__cobra_collisions[vidx]
            cobra_assignments = self.__cobra_assignments[vidx]
            target_assignments = self.__target_assignments[vidx]

            # Get the target ids for each colliding cobra
            tidx = cobra_assignments[collisions]

            # Look up the target classes and non-observation costs for non-broken cobras
            mask = tidx != -1
            target_class = self.__target_cache.target_class[tidx[mask]]
            cost = np.full_like(tidx, -1)
            for i, ix in enumerate(zip(*np.where(mask))):
                if target_class[i] in self.__netflow_options.target_classes:
                    cost[ix] = self.__netflow_options.target_classes[target_class[i]].non_observation_cost
                else:
                    # This is a target with missing priority
                    cost[ix] = 0
            
            # Unassign the target with the lower cost, except when one of the cobras is broken
            # because in that case the good cobra should be unassigned as it would collide with
            # the broken cobra
            cost[cost == -1] = np.max(cost) + 1
            cidx = collisions[np.argmin(cost, axis=0), np.arange(collisions.shape[1])]

            # First remove the target assignments
            for ci in cidx:
                target_assignments.pop(cobra_assignments[ci], None)

            # Then unassign the cobras
            cobra_assignments[cidx] = -1

            if collisions.size > 0:
                logger.info(f'Unassigned cobras {cidx + 1} for visit {vidx} due to collisions.')
        
    def __get_fiber_status(self):
        """
        Return the fiber status indexes by cidx, as required by PfsDesign.
        """

        # 1-based cobra IDs, should have the size of 2394
        cobraid = np.arange(self.__bench.cobras.nCobras, dtype=int) + 1

        # 0-based indices of science fibers, in order of cobraId
        sci_fiberidx = self.__fiber_map.cobraIdToFiberId(cobraid) - 1

        fiber_status = np.full(self.__bench.cobras.nCobras, FiberStatus.GOOD, dtype=np.int32)

        fiber_broken_mask = (self.__calib_model.status & self.__calib_model.FIBER_BROKEN_MASK) != 0
        cobra_ok_mask = (self.__calib_model.status & self.__calib_model.COBRA_OK_MASK) != 0
        cobra_broken_mask = ~cobra_ok_mask & ~fiber_broken_mask
        fiber_blocked_mask = np.array(self.__blocked_fibers.loc[self.__fiber_map.fiberId[sci_fiberidx]].status)

        fiber_status[fiber_broken_mask] = FiberStatus.BROKENFIBER
        fiber_status[cobra_broken_mask] = FiberStatus.BROKENCOBRA
        fiber_status[fiber_blocked_mask] = FiberStatus.BLOCKED

        return fiber_status

    def get_fiber_assignments(self,  
                              include_target_columns=False,
                              include_unassigned_fibers=False,
                              include_engineering_fibers=False):
        
        """
        Return a data frame of the target assignments including positions, visitid, fiberid
        and other information.

        When `all_columns` is set to True, all columns from the target catalog are included in the
        output.

        Note that column names ending with `_idx` are 0-based indices into netflow data structured,
        while the columns `fiberid` and `cobraid` are 1-based indices, as defined by PFS.

        This is very similar what happens in pfs.utils.pfsDesignUtils.makePfsDesign but we do it for
        all visits at once and work on a DataFrame instead of PfsDesign objects and don't assume a
        specific ordering of the fibers or the cobras in the DataFrame.

        Parameters:
        -----------
        include_target_columns : bool
            Include all columns from the input target catalog in the output.
        include_unassigned_fibers : bool
            Include unassigned fibers in the output.
        include_engineering_fibers : bool
            Include engineering fibers in the output.
        include_empty_fibers : bool
            Include empty fibers in the output.
        """

        # TODO: this function uses a lot of PFS specific logic and calibration products
        #       consider moving it elsewhere, maybe to the Design class which currently only
        #       converts between data frames and pfsDesign objects

        fm = self.__fiber_map

        # There are 2394 cobras in total
        # There are 2604 fibers, 2394 assigned to cobras, 64 engineering fibers and 146 empty fibers
        # The design files contain 2458 rows, for the cobras plus the 64 engineering fibers

        # Internally, cidx is a 0-based index of the cobras, PFS cobraids are 1-based
        # Map the cobraIds to the corresponding fiberIds
        # These include all cobras that were part of the netflow problem

        # TODO: make sure this is updated if we leave out cobras from the netflow problem
        #       because they're broken or else

        # Map cidx (zero-based) to various science fiber identitfiers
        # Should have the size of 2394
        cobraid = np.arange(self.__bench.cobras.nCobras, dtype=int) + 1
        sci_fiberidx = fm.cobraIdToFiberId(cobraid) - 1                     # Indices of science fibers, in order of cobraId
        eng_fiberidx = np.where(fm.scienceFiberId == fm.ENGINEERING)[0]     # Indices of engineering fibers, arbitrary order

        # TODO: consider removing broken fibers from the netflow problem at the beginning
        #       instead of updating the status now
        fiber_status = self.__get_fiber_status()

        assignments : pd.DataFrame = None
        
        for vidx, visit in enumerate(self.__visits):
            # Unassigned cobras
            cidx = np.where(self.__cobra_assignments[vidx] == -1)[0]
            stage = visit.pointing.stage if visit.pointing.stage is not None else -1

            if include_unassigned_fibers:
                unassigned = pd.DataFrame({ 'fiberid': fm.fiberId[sci_fiberidx][cidx] })
                # pd_append_column(unassigned, 'key', None, 'string')
                # pd_append_column(unassigned, 'targetid', -1, np.int64)
                pd_append_column(unassigned, 'stage', stage, np.int32)
                pd_append_column(unassigned, 'pointing_idx', visit.pointing_idx, np.int32)
                pd_append_column(unassigned, 'visit_idx', vidx, np.int32)
                pd_append_column(unassigned, 'target_idx', -1, np.int32)
                pd_append_column(unassigned, 'cobraid', cobraid[cidx], np.int32)
                pd_append_column(unassigned, 'sciencefiberid', fm.scienceFiberId[sci_fiberidx][cidx], np.int32)
                pd_append_column(unassigned, 'fieldid', fm.fieldId[sci_fiberidx][cidx], np.int32)
                pd_append_column(unassigned, 'fiberholeid', fm.fiberHoleId[sci_fiberidx][cidx], np.int32)
                pd_append_column(unassigned, 'spectrographid', fm.spectrographId[sci_fiberidx][cidx], np.int32)
                pd_append_column(unassigned, 'fp_x', self.__bench.cobras.home0[cidx].real, np.float64)
                pd_append_column(unassigned, 'fp_y', self.__bench.cobras.home0[cidx].imag, np.float64)
                pd_append_column(unassigned, 'target_type', 'na', 'string')
                pd_append_column(unassigned, 'fiber_status', fiber_status[cidx], np.int32)
                pd_append_column(unassigned, 'center_dist', np.nan, np.float64)
                pd_append_column(unassigned, 'black_dot_dist', np.nan, np.float64)
                
                if assignments is None:
                    assignments = unassigned
                else:
                    assignments = pd.concat([ assignments, unassigned ])

            if include_engineering_fibers:
                engineering = pd.DataFrame({ 'fiberid': fm.fiberId[eng_fiberidx] })
                # pd_append_column(engineering, 'key', None, 'string')
                # pd_append_column(engineering, 'targetid', -1, np.int64)
                pd_append_column(engineering, 'stage', stage, np.int32)
                pd_append_column(engineering, 'pointing_idx', visit.pointing_idx, np.int32)
                pd_append_column(engineering, 'visit_idx', vidx, np.int32)
                pd_append_column(engineering, 'target_idx', -1, np.int32)
                pd_append_column(engineering, 'cobraid', -1, np.int32)
                pd_append_column(engineering, 'sciencefiberid', fm.scienceFiberId[eng_fiberidx], np.int32)
                pd_append_column(engineering, 'fieldid', fm.fieldId[eng_fiberidx], np.int32)
                pd_append_column(engineering, 'fiberholeid', fm.fiberHoleId[eng_fiberidx], np.int32)
                pd_append_column(engineering, 'spectrographid', fm.spectrographId[eng_fiberidx], np.int32)
                pd_append_column(engineering, 'fp_x', np.nan, np.float64)
                pd_append_column(engineering, 'fp_y', np.nan, np.float64)
                pd_append_column(engineering, 'target_type', 'eng', 'string')
                pd_append_column(engineering, 'fiber_status', FiberStatus.GOOD, np.int32)
                pd_append_column(engineering, 'center_dist', np.nan, np.float64)
                pd_append_column(engineering, 'black_dot_dist', np.nan, np.float64)

                if assignments is None:
                    assignments = engineering
                else:
                    assignments = pd.concat([ assignments, engineering ])

            # Assigned targets
            tidx = np.array([ k for k in self.__target_assignments[vidx] ], dtype=int)      # target index
            cidx = np.array([ self.__target_assignments[vidx][ti] for ti in tidx ])         # cobra index, 0-based
            fpidx = self.__cache_to_fp_pos_map[visit.pointing_idx][tidx]
            fp_pos = self.__target_fp_pos[visit.pointing_idx][fpidx]

            # Calculate the distance from the cobra center and the closes black dots for each target
            center_dist = np.full_like(cidx, np.nan, dtype=np.float64)
            black_dot_dist = np.full_like(cidx, np.nan, dtype=np.float64)
            for i, ci in enumerate(cidx):
                center_dist[i] = self.__get_cobra_center_distance(ci, fp_pos[i])
                black_dot_dist[i] = self.__get_closest_black_dot_distance(ci, fp_pos[i])

            targets = pd.DataFrame({ 'fiberid': fm.fiberId[sci_fiberidx][cidx] })
            # pd_append_column(targets, 'key', self.__target_cache.key[tidx], 'string')
            # pd_append_column(targets, 'targetid', self.__target_cache.id[tidx], np.int64)
            pd_append_column(targets, 'stage', stage, np.int32)
            pd_append_column(targets, 'pointing_idx', visit.pointing_idx, np.int32)
            pd_append_column(targets, 'visit_idx', vidx, np.int32)
            pd_append_column(targets, 'target_idx', self.__target_cache.target_idx[tidx], np.int32)
            pd_append_column(targets, 'cobraid', cobraid[cidx], np.int32)
            pd_append_column(targets, 'sciencefiberid', fm.scienceFiberId[sci_fiberidx][cidx], np.int32)
            pd_append_column(targets, 'fieldid', fm.fieldId[sci_fiberidx][cidx], np.int32)
            pd_append_column(targets, 'fiberholeid', fm.fiberHoleId[sci_fiberidx][cidx], np.int32)
            pd_append_column(targets, 'spectrographid', fm.spectrographId[sci_fiberidx][cidx], np.int32)
            pd_append_column(targets, 'fp_x', fp_pos.real, np.float64)
            pd_append_column(targets, 'fp_y', fp_pos.imag, np.float64)
            pd_append_column(targets, 'target_type', self.__target_cache.prefix[tidx], 'string')
            pd_append_column(targets, 'fiber_status', fiber_status[cidx], np.int32)
            pd_append_column(targets, 'center_dist', center_dist, np.float64)
            pd_append_column(targets, 'black_dot_dist', black_dot_dist, np.float64)

            if assignments is None:
                assignments = targets
            else:
                assignments = pd.concat([ assignments, targets ])

        # Map target prefix to PFS target type
        # Map internal prefixes to PFS fiber status
        target_type_map = {
            'na': TargetType.UNASSIGNED,
            'sky': TargetType.SKY,
            'cal': TargetType.FLUXSTD,
            'sci': TargetType.SCIENCE,
            'eng': TargetType.ENGINEERING,
        }

        assignments['target_type'] = assignments['target_type'].map(target_type_map).astype(np.int32)

        # Include all columns from the target list data frame, if requested
        # Convert integer columns to nullable to avoid float conversion in join
        if include_target_columns:  
            assignments = assignments.join(pd_to_nullable(self.__targets), on='target_idx', how='left')

        # Make sure float columns contain NaN instead of None
        pd_null_to_nan(assignments, in_place=True)

        return assignments
    
    def get_fiber_assignments_masks(self):
        """
        Returns a list of masks that indexes the assigned targets within the target list.
        """

        masks = []
        fiberids = []
        for i, p in enumerate(self.__visits):
            tidx = self.__cache_to_target_map[np.array(list(self.__target_assignments[i].keys()))]
            mask = np.full(len(self.__targets), False, dtype=bool)
            mask[tidx] = True
            masks.append(mask)
            fiberids.append(np.array(self.__target_assignments[i].values()))

        return masks, fiberids
    
    def get_target_assignment_summary(self):
        """
        Return a data frame of the input targets with the number of visits each target was
        scheduled for.

        This output is useful for identifying targets that were not assigned to any visit or
        targets that were assigned to more than the required number of visits.
        """

        all_assignments, all_fiberids = self.get_fiber_assignments_masks()
        num_assignments = np.sum(np.stack(all_assignments, axis=-1), axis=-1)

        targets = self.__targets.copy()
        targets['num_visits'] = num_assignments

        return targets

    def __get_idcol(self, catalog):
        for c in ['objid', 'skyid', 'calid']:
            if c in catalog.data.columns:
                return c
        
        raise NetflowException()