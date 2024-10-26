import os
import re
from typing import Callable
from collections import defaultdict
from collections.abc import Iterable
import numpy as np
import pandas as pd
from types import SimpleNamespace  

from pfs.datamodel import TargetType, FiberStatus
from .setup_logger import logger
from .util import *
from ..util.args import *
from ..util.config import *
from ..projection import Pointing
from .visit import Visit
from .gurobiproblem import GurobiProblem

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
        * elbow_collisions : bool
            If True, check for collisions between the fiber elbows, otherwise check for collisions
            between the fiber ends.
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
        * ignore_forbidden_pairs
        * ignore_forbidden_singles
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
    * ID          - must be unique across all targets
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
    
    The problem can be solved with the `solve` function which executes the solver and calls
    `__extract_assignments` to extract the results. The results are stored in the variables
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
        self.__netflow_options = camel_to_snake(netflow_options)
        self.__debug_options = camel_to_snake(debug_options)

        # Internally used variables

        self.__bench = instrument.bench             # Bench object
        self.__fiber_map = instrument.fiber_map     # Grand fiber map

        self.__problem_type = GurobiProblem
        self.__problem = None                       # ILP problem, already wrapped

        self.__name_counter = None

        self.__variables = None
        self.__constraints = None

        self.__target_classes = None
        self.__forbidden_targets = None             # List of forbidden individual target, identified by target_idx
        self.__forbidden_pairs = None               # List of forbidden target pairs, identified by target_idx
        self.__black_dots = None
        self.__cobra_groups = None
        self.__time_budgets = None

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
        """
        Return an option from the netflow_options dictionary, if it exists,
        otherwise return the default value.
        """

        if self.__netflow_options is not None and key in self.__netflow_options:
            return self.__netflow_options[key]
        else:
            return default
        
    def __get_debug_option(self, key, default=None):
        """
        Return an option from the debug_options dictionary, if it exists,
        otherwise return the default value
        """

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
        self.__extract_assignments()

    def save_solution(self, filename=None):
        """Save the LP solution to a file."""

        fn = self.__get_solution_filename(filename)
        self.__problem.write_solution(fn)

    #endregion
    #region Configuration
        
    def __get_target_class_config(self):
        target_classes = {}

        for name, options in self.__get_netflow_option('target_classes', {}).items():
            
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
    
    def __get_forbidden_targets_config(self):
        """
        Look up the target index based on the target id of forbidden targets.
        """

        forbidden_targets = self.__get_netflow_option('forbidden_targets', None)
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

        forbidden_pairs = self.__get_netflow_option('forbidden_pairs', None)
        fpp = []
        wrong_id_count = 0
        if forbidden_pairs is not None:
            for i, pair in enumerate(forbidden_pairs):
                if len(pair) != 2:
                    raise ValueError(f"Found an incorrect number of target ids in forbidden pair list at index {i}.")
            
                tidx_list = [ self.__target_to_cache_map[self.__targets.index.get_loc(p)] for p in pair ]
                if -1 not in tidx_list:
                    fpp.append(tidx_list)
                else:
                    wrong_id_count += 1

        if wrong_id_count:
            logger.warning(f"Found {wrong_id_count} forbidden targets pairs that are outside the pointings.")
                
        return fpp

    def __get_black_dot_config(self):
        black_dot_penalty = self.__get_netflow_option('black_dot_penalty', None)
        
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
        
        for name, options in self.__get_netflow_option('cobra_groups', {}).items():
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

        for name, options in self.__get_netflow_option('time_budgets', {}).items():
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

    def __calculate_fp_pos(self, pointing, ra, dec, pmra=None, pmdec=None, parallax=None, rv=None):

        epoch = self.__get_netflow_option('epoch', 2016.0)
        ignore_proper_motion = self.__get_netflow_option('ignore_proper_motion', False)

        # The proper motion of the targets used for sky_pfi transformation.
        # The unit is mas/yr, the shape is (2, N)
        if ignore_proper_motion:
            pmra = pmdec = None

        # The parallax of the coordinates used for sky_pfi transformation.
        # The unit is mas, the shape is (1, N)
        if ignore_proper_motion:
            parallax = None
            pmra = pmdec = None
            rv = None

        fp_pos = self.__instrument.radec_to_fp_pos(pointing, ra, dec,
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

        res = []
        for cidx in range(len(self.__bench.cobras.centers)):
            nb = self.__bench.getCobraNeighbors(cidx)
            res.append(self.__bench.blackDots.centers[[cidx] + list(nb)])

        return res
    
    def __get_visibility_and_elbow(self, fp_pos, fpidx_to_tidx_map):
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
            for i, fpidx in enumerate(targets[cidx, :]):
                if fpidx >= 0:
                    # Map focal plane index to the target cache index
                    tidx = fpidx_to_tidx_map[fpidx]
                    targets_cobras[tidx].append((cidx, elbows[cidx, i]))
                    cobras_targets[cidx].append((tidx, elbows[cidx, i]))

        return targets_cobras, cobras_targets
    
    def __get_colliding_pairs(self, pidx, collision_distance):
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
        set : pairs of target indices that would cause fiber top collisions.
        """

        fp_pos = self.__target_fp_pos[pidx]
        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]
        visibility = self.__visibility[pidx]

        # Collect targets associated with each cobra, for each visit
        fp_pos = np.array(fp_pos)
        ivis = defaultdict(list)
        for tidx, cidx_elbow in visibility.targets_cobras.items():
            for (cidx, _) in cidx_elbow:
                ivis[cidx].append(tidx)

        pairs = set()
        for cidx, tidx1 in ivis.items():
            # Determine target indices visible by this cobra and its neighbors
            nb = self.__bench.getCobraNeighbors(cidx)
            tidx2 = np.concatenate([ivis[j] for j in nb if j in ivis])
            tidx2 = np.concatenate((tidx1, tidx2))
            tidx2 = np.unique(tidx2).astype(int)
            d = np.abs(np.subtract.outer(fp_pos[tidx_to_fpidx_map[tidx1]], fp_pos[tidx_to_fpidx_map[tidx2]]))
            
            for m, n in zip(*np.where(d  < collision_distance)):
                if tidx1[m] < tidx2[n]:               # Only store pairs once
                    pairs.add((tidx1[m], tidx2[n]))

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

        visibility = self.__visibility[pidx]
        fp_pos = self.__target_fp_pos[pidx]
        tidx_to_fpidx_map = self.__cache_to_fp_pos_map[pidx]

        res = defaultdict(list)
        for cidx, tidx_elbow in visibility.cobras_targets.items():
            # Determine target indices visible by neighbors of this cobra
            
            # List of neighboring cobra_index
            nb = self.__bench.getCobraNeighbors(cidx)       

            # All targets and elbow positions visible by neighboring cobras
            # tmp = [ vis_cobras_targets[j] for j in nb if j in vis_cobras_targets ]
            tmp = [ ti for j in nb if j in visibility.cobras_targets for ti, _ in visibility.cobras_targets[j] ]

            if len(tmp) > 0:
                # Unique list of target indices visible by neighboring cobras
                i2 = np.unique(tmp)

                # For each target visible by this cobra and the corresponding elbow
                # position, find all targets which are too close to the "upper arm"
                # of the cobra
                for tidx, elbowpos in tidx_elbow:
                    ebp = np.full(len(i2), elbowpos)                            # elbow position of the cobra
                    tp = np.full(len(i2), fp_pos[tidx_to_fpidx_map[tidx]])      # target position
                    ti2 = fp_pos[tidx_to_fpidx_map[i2]]                         # target positions of the neighbors
                    d = self.__bench.distancesToLineSegments(ti2, tp, ebp)
                    res[(cidx, tidx)] += list(i2[d < collision_distance])

        return res

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
        
    def __append_targets(self, catalog, id_column, prefix, exp_time=None, priority=None, penalty=None, mask=None, filter=None, selection=None):
        
        # Make sure that we convert to the right data type everywhere
                       
        df = catalog.get_data(mask=mask, filter=filter, selection=selection)

        # Create the formatted dataset from the input catalog
        data = {
            'id': df[id_column].astype(np.int64).reset_index(drop=True),
            'RA': df['RA'].astype(np.float64).reset_index(drop=True),
            'Dec': df['Dec'].astype(np.float64).reset_index(drop=True)
        }
        targets = pd.DataFrame(data)

        # Add proper motion and parallax, if available
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
        
        # This is a per-object penalty for observing calibration targets
        if penalty is not None:
            pd_append_column(targets, 'penalty', penalty, np.int32)
        elif 'penalty' in df.columns:
            pd_append_column(targets, 'penalty', df['penalty'], np.int32)
        else:
            pd_append_column(targets, 'penalty', 0, np.int32)

        pd_append_column(targets, 'prefix', prefix, str)

        if prefix == 'sci':
            if exp_time is not None:
                pd_append_column(targets, 'exp_time', exp_time, np.float64)
            else:
                pd_append_column(targets, 'exp_time', df['exp_time'], np.float64)

            if priority is not None:
                pd_append_column(targets, 'priority', priority, np.int32)
            else:
                pd_append_column(targets, 'priority', df['priority'], np.int32)

            targets['class'] = targets[['prefix', 'priority']].apply(lambda r: f"{r['prefix']}_P{r['priority']}", axis=1)
        else:
            pd_append_column(targets, 'class', prefix, str)

            # Calibration targets have no prescribed exposure time and priority
            pd_append_column(targets, 'exp_time', np.nan, np.float64)
            pd_append_column(targets, 'priority', -1, np.int32)

        # Append to the existing list
        if self.__targets is None:
            self.__targets = targets
        else:
            self.__targets = pd.concat([self.__targets, targets])
            
    def append_science_targets(self, catalog, exp_time=None, priority=None, mask=None, filter=None, selection=None):
        """Add science targets"""

        self.__append_targets(catalog, 'objid', 'sci', exp_time=exp_time, priority=priority, mask=mask, filter=filter, selection=selection)

    def append_sky_targets(self, sky, mask=None, filter=None, selection=None):
        """Add sky positions"""

        self.__append_targets(sky, 'skyid', prefix='sky', mask=mask, filter=filter, selection=selection)

    def append_fluxstd_targets(self, fluxstd, mask=None, filter=None, selection=None):
        """Add flux standard positions"""

        self.__append_targets(fluxstd, 'objid', prefix='cal', mask=mask, filter=filter, selection=selection)

    def build(self):
        """Construct the ILP problem"""
        # Load configuration
        logger.info("Processing netflow configuration")

        # Run a few sanity checks
        self.__check_pointing_visibility()

        # Optimize data access
        self.__calculate_exp_time()
        self.__calculate_target_visits()
        self.__cache_targets()

        # Target classes
        self.__target_classes = self.__get_target_class_config()

        # Forbidden targets and target pairs
        self.__forbidden_targets = self.__get_forbidden_targets_config()
        self.__forbidden_pairs = self.__get_forbidden_pairs_config()

        # Cobras positioned too close to a black dot can get a penalty
        self.__black_dots = self.__get_black_dot_config()

        # Cobra groups are defined to set a minimum number of calibration targets in each
        self.__cobra_groups = self.__get_cobra_groups_config()

        # Science program time budgets
        self.__time_budgets = self.__get_time_budget_config()

        # Build the problem
        self.__create_visits()
        self.__calculate_target_fp_pos()

        # Run a few sanity checks
        self.__check_target_fp_pos()

        self.__build_ilp_problem()

    def __check_pointing_visibility(self):
        """
        Verify that all pointings are visible in the sky at obs_time.
        """

        for p in self.__pointings:
            # Convert pointing center into azimuth and elevation (+ instrument rotator angle)
            alt, az, inr = self.__instrument.radec_to_altaz(p.ra, p.dec, p.posang, p.obs_time)
            assert az > 0, f"Pointing is below the horizon."
                    
    def __calculate_exp_time(self):
        """
        Calculate the exposure time for each object.
        """

        # Since targeting is done in fix quanta in time, make sure all pointings and visits
        # have the same exposure time

        self.__visit_exp_time = None
        for p in self.__pointings:
            if self.__visit_exp_time is None:
                self.__visit_exp_time = p.exp_time
            elif self.__visit_exp_time != p.exp_time:
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
        self.__targets.loc[sci, 'req_visits'] = np.int32(np.ceil(self.__targets['exp_time'][sci] / self.__visit_exp_time))

        max_req_visits = self.__targets['req_visits'].max()
        for p in self.__pointings:
            if p.nvisits < max_req_visits:
                raise RuntimeError('Some science targets require more visits than provided.')

    def __cache_targets(self):
        """Extract the contents of the Pandas DataFrame for faster indexed access."""

        # Sort targets by id, if not already
        if 'id' in self.__targets.columns:
            self.__targets = self.__targets.set_index('id').sort_index()

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
                id = np.array(self.__targets.index[mask_any].astype(np.int64)),
                ra = np.array(self.__targets['RA'][mask_any].astype(np.float64)),
                dec = np.array(self.__targets['Dec'][mask_any].astype(np.float64)),
                pm = np.array(self.__targets['pm'][mask_any].astype(np.float64)),
                pmra = np.array(self.__targets['pmra'][mask_any].astype(np.float64)),
                pmdec = np.array(self.__targets['pmdec'][mask_any].astype(np.float64)),
                parallax = np.array(self.__targets['parallax'][mask_any].astype(np.float64)),
                target_class = np.array(self.__targets['class'][mask_any].astype(str)),
                prefix = np.array(self.__targets['prefix'][mask_any].astype(str)),
                req_visits = np.array(self.__targets['req_visits'][mask_any].astype(np.int32)),
                done_visits = np.array(self.__targets['done_visits'][mask_any].astype(np.int32)),
                penalty = np.array(self.__targets['penalty'][mask_any].astype(np.int32)),
            )

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
                                             parallax=self.__target_cache.parallax[mask])
            
            self.__target_fp_pos.append(fp_pos)

    def __check_target_fp_pos(self):
        # Verify if each pointing contains a reasonable number of targets
        # A target is accessible if it is within the 400mm radius of the PFI

        for fp in self.__target_fp_pos:
            n = np.sum(np.abs(fp) < 400)
            assert n > 100, f"Pointing contains only {n} targets within the PFI radius."

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
        if isinstance(constraint, tuple):
            self.__problem.add_linear_constraint(name, *constraint)
        else:
            self.__problem.add_constraint(name, constraint)

    def __build_ilp_problem(self):
        """
        Construct the ILP problem by defining the variables and constraints.
        """

        use_named_variables = self.__get_netflow_option('use_named_variables', False)
        if use_named_variables:
            self.__make_name = self.__make_name_full
        else:
            self.__name_counter = 0
            self.__make_name = self.__make_name_short

        # Number of visits
        nvisits = len(self.__visits)

        self.__init_problem()
        self.__init_variables()
        self.__init_constraints()

        logger.info("Calculating target visibilities")
        self.__calculate_visibilities()

        logger.info("Calculating cobra collisions")
        self.__calculate_collisions()

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

            # > Tv_o_coll_{?}_{?}_{visit_idx} (endpoint collisions)
            # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx} (elbow collisions)
            logger.debug("Creating cobra collision constraints.")
            self.__create_cobra_collision_constraints(visit)

            logger.debug("Adding cobra non-allocation cost terms.")
            self.__add_cobra_non_allocation_cost(visit)

            # Add constraints for forbidden targets and forbidden pairs
            ignore_forbidden_targets = self.__get_debug_option('ignore_forbidden_targets', False)
            if not ignore_forbidden_targets:
                if self.__forbidden_targets is not None and len(self.__forbidden_targets) > 0:
                    # > Tv_o_forb_{tidx}_{tidx}_{visit_idx}
                    logger.debug("Adding forbidden target constraints")
                    self.__create_forbidden_target_constraints(visit)
            else:
                 logger.debug("Ignored forbidden target constraints")

            ignore_forbidden_pairs = self.__get_debug_option('ignore_forbidden_pairs', False)
            if not ignore_forbidden_pairs:
                if self.__forbidden_pairs is not None and len(self.__forbidden_pairs) > 0:
                    # > Tv_o_forb_{tidx1}_{tidx2}_{visit_idx}
                    logger.debug("Adding forbidden pair constraints")
                    self.__create_forbidden_pair_constraints(visit)
            else:
                logger.debug("Ignroed forbidden pair constraints")
        
        logger.info("Adding constraints")

        # TODO: add log message in function beyond this point

        # The total number of science targets per target class must be balanced and
        # the minimum and maximum number of targets to be observed within
        # a science target class can be optionally constrained.
        # > STC_o_sum_{target_class}
        # > STC_o_min_{target_class}
        # > STC_o_max_{target_class}
        self.__create_science_target_class_constraints()

        # Every calibration target class must be observed a minimum number of times
        # every visit
        # > CTCv_o_sum_{target_class}_{visit_idx}
        # > CTCv_o_min_{target_class}_{visit_idx}
        # > CTCv_o_max_{target_class}_{visit_idx}
        self.__create_calibration_target_class_constraints()

        # Inflow and outflow at every T node must be balanced
        # > T_i_T_o_sum_{target_id}
        self.__create_science_target_constraints()

        # Inflow and outflow at every Tv node must be balanced
        # > Tv_i_Tv_o_sum_{target_id}_{visit_idx}
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

    def __calculate_visibilities(self):
        self.__visibility = []
        for pidx, pointing in enumerate(self.__pointings):
            logger.info(f'Calculating visibility for pointing {pidx + 1}.')

            # Get the target visibility and elbow positions for each target roughly within the pointing
            fp_pos = self.__target_fp_pos[pidx]
            fpidx_to_tidx_map = self.__fp_pos_to_cache_map[pidx]

            vis = SimpleNamespace()
            vis.targets_cobras, vis.cobras_targets = self.__get_visibility_and_elbow(fp_pos, fpidx_to_tidx_map)
            self.__visibility.append(vis)

            # targets_cobras keys   -> index into the target_cache
            # targets_cobras values -> (cidx, elbow position)
            
            # cobras_targets keys   -> cobra index
            # cobras_targets values -> (tidx, elbow position), where tidx is the index into the target_cache

    def __calculate_collisions(self):
        collision_distance = self.__get_netflow_option('collision_distance', 0.0)

        self.__collisions = []
        for pidx, pointing in enumerate(self.__pointings):
            logger.info(f'Calculating collisions for pointing {pidx + 1}.')

            # Get the target visibility and elbow positions for each target roughly within the pointing
            col = SimpleNamespace()
            col.pairs = self.__get_colliding_pairs(pidx, collision_distance)
            col.elbows = self.__get_colliding_elbows(pidx, collision_distance)
            self.__collisions.append(col)

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
        """
        Create the variables corresponding to the STC -> T and STC -> STC_sink edges.
        These edges are created regardless of the visibility of the targets during the visits.
        """

        # TODO: replace this with the logic to implement step-by-step targeting
        force_already_observed = self.__get_netflow_option('force_already_observed', False)

        # Select only science targets and create the variables in batch mode
        mask = self.__target_cache.prefix == 'sci'
        tidx = np.where(mask)[0]
        target_class = self.__target_cache.target_class[mask]
        target_penalty = self.__target_cache.penalty[mask]

        # Science Target class node to target node: STC -> T
        self.__create_STC_T(tidx, target_class, target_penalty)

        # Science Target node to sink: T -> T_sink
        # > T_sink[target_idx]
        self.__create_T_sink(tidx, target_class)

    def __create_STC_T(self, tidx, target_class, target_penalty):
        vars =  self.__add_variable_array('STC_T', tidx, 0, 1)
        for i in range(len(tidx)):
            f = vars[tidx[i]]
            self.__variables.T_i[tidx[i]] = f
            self.__variables.STC_o[target_class[i]].append(f)
            if target_penalty[i] != 0:
                self.__add_cost(f * target_penalty[i])
        
    def __create_T_sink(self, tidx, target_class):
        """
        Create the sink variable which is responsible for draining the flow from the target nodes
        which get less visits than required
        """
        max_visits = len(self.__visits)
        T_sink =  self.__add_variable_array('T_sink', tidx, 0, max_visits)

        for i in range(len(tidx)):
            self.__variables.T_sink[tidx[i]] = T_sink[tidx[i]]
            self.__add_cost(T_sink[tidx[i]] * self.__target_classes[target_class[i]].partial_observation_cost)

    def __create_calib_target_class_variables(self, visit):
        # Calibration target class visit outflows, key: (target_class, visit_idx)
        for target_class, options in self.__target_classes.items():
            if options.prefix in ['sky', 'cal']:
                f = self.__add_variable(self.__make_name("CTCv_sink", target_class, visit.visit_idx), 0, None)
                self.__variables.CTCv_sink[(target_class, visit.visit_idx)] = f
                self.__add_cost(f * options.non_observation_cost)
    
    def __create_visit_variables(self, visit):
        """
        Create the variables for the visible variables for a visit.
        """

        visibility = self.__visibility[visit.pointing_idx]
        tidx = np.array(list(visibility.targets_cobras.keys()), dtype=int)

        # Create science targets T_Tv edges in batch
        # > T_Tv_{visit_idx}[target_idx]
        sci_mask = self.__target_cache.prefix[tidx] == 'sci'
        self.__create_T_Tv(visit, tidx[sci_mask])

        # Create calibration targets CTCv_Tv edges
        # > CTCv_Tv_{visit_idx}[target_idx]
        cal_mask = (self.__target_cache.prefix[tidx] == 'sky') | (self.__target_cache.prefix[tidx] == 'cal')
        self.__create_CTCv_Tv(visit, tidx[cal_mask])

        # For each Cv, generate the incoming Tv_Cv edges
        # These differ in number but can be vectorized for each cobra
        # > Tv_Cv_{visit_idx}_{cobra_idx}[target_idx]
        for cidx, tidx_elbow in visibility.cobras_targets.items():
            tidx = np.array([ ti for ti, _ in tidx_elbow ], dtype=int)
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

        # Cost of moving the cobra away from the center
        cobra_move_cost = self.__get_netflow_option('cobra_move_cost', None)
        if cobra_move_cost is not None:
            dist = np.abs(self.__bench.cobras.centers[cidx] - 
                          self.__target_fp_pos[visit.pointing_idx][tidx_to_fpidx_map[tidx]])
            cost += cobra_move_cost(dist)
        
        # Cost of closest black dots for each cobra
        if self.__black_dots is not None:
            # Distance of all targets within the patrol radius to the nearby black dots
            dist = np.abs(self.__black_dots.black_dot_list[cidx] - 
                          self.__target_fp_pos[visit.pointing_idx][tidx_to_fpidx_map[tidx]][:, None])
            # Take minimum distance from each target
            cost += self.__black_dots.black_dot_penalty(np.min(dist, axis=-1))

        # Create LP variables
        name = self.__make_name("Tv_Cv", visit.visit_idx, cidx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        # For each accessible target, save the Tv -> Cv variable to the respective lists
        # TODO: can these be created in batch mode?
        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]

            # TODO: consider creating a separate list for each visit
            self.__variables.Cv_i[(cidx, visit.visit_idx)].append(f)
            self.__variables.Tv_o[(ti, visit.visit_idx)].append((f, cidx))

            # Save the variable to the list of each cobra group to which it's relevant
            for cg_name, options in self.__cobra_groups.items():
                if target_class in options.target_classes:
                    self.__variables.CG_i[(cg_name, visit.visit_idx, options.groups[cidx])].append(f)

    #endregion
    #region Special cost term
            
    def __add_cobra_non_allocation_cost(self, visit):
        # If requested, penalize non-allocated fibers
        # Sum up the all Cv_i edges for the current visit and penalize its difference from
        # the total number of cobras
        cobra_non_allocation_cost = self.__get_netflow_option('fiber_non_allocation_cost', 0)
        if cobra_non_allocation_cost != 0:
            # TODO: consider storing Cv_i organized by visit instead of cobra as well
            relevant_vars = [ var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == visit.visit_idx ]
            relevant_vars = [ item for sublist in relevant_vars for item in sublist ]
            self.__add_cost(cobra_non_allocation_cost *
                            (self.__bench.cobras.nCobras - self.__problem.sum(relevant_vars)))
            
    #endregion
    #region Constraints

    def __create_cobra_collision_constraints(self, visit):
        ignore_endpoint_collisions = self.__get_debug_option('ignore_endpoint_collisions', False)
        ignore_elbow_collisions = self.__get_debug_option('ignore_elbow_collisions', False)
        
        # Add constraints 
        logger.debug(f"Adding constraints for visit {visit.visit_idx}")

        # Avoid endpoint or elbow collisions
        collision_distance = self.__get_netflow_option('collision_distance', 0.0)
        elbow_collisions = self.__get_netflow_option('elbow_collisions', False)
        
        if collision_distance > 0.0:
            # if not elbow_collisions:
                if not ignore_endpoint_collisions:
                    # > Tv_o_coll_{?}_{?}_{visit_idx}
                    logger.debug("Adding endpoint collision constraints")
                    self.__create_endpoint_collision_constraints(visit)
                else:
                    logger.debug("Ignoring endpoint collision constraints")
            # else:
                if not ignore_elbow_collisions:
                    # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx}
                    logger.debug("Adding elbow collision constraints")
                    self.__create_elbow_collision_constraints(visit)
                else:
                    logger.debug("Ignoring elbow collision constraints")
                
    def __create_endpoint_collision_constraints(self, visit):
        vidx = visit.visit_idx
        collisions = self.__collisions[visit.pointing_idx]

        # TODO: delete
        # Collect all targets that are relevant for this visit
        # TODO: this is not the most optimal because have to iterate over all visits
        # tidx = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == vidx)

        # TODO: delete
        # keys = self.__variables.Tv_o.keys()
        # keys = set(key[0] for key in keys if key[1] == vidx)
        
        for p0, p1 in collisions.pairs:
            # TODO: is this check necessary?
            # if p0 in tidx and p1 in tidx:

            vars = [ v for v, cidx in self.__variables.Tv_o[(p0, vidx)] ] + \
                   [ v for v, cidx in self.__variables.Tv_o[(p1, vidx)] ]
            
            name = self.__make_name("Tv_o_coll", p0, p1, vidx)
            # constr = self.__problem.sum(vars) <= 1
            constr = ([1] * len(vars), vars, '<=', 1)
            self.__constraints.Tv_o_coll[(p0, p1, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_elbow_collision_constraints(self, visit):
        vidx = visit.visit_idx
        collisions = self.__collisions[visit.pointing_idx]

        for (cidx1, tidx1), tidx2_list in collisions.elbows.items():

            # TODO: inefficient loop
            for f, cidx2 in self.__variables.Tv_o[(tidx1, vidx)]:
                if cidx2 == cidx1:
                    var0 = f

            for tidx2 in tidx2_list:
                # TODO: is this test necessary
                # if True:  # idx2 != tidx1:
                
                vars = [ var0 ]
                vars += [ f for f, cidx2 in self.__variables.Tv_o[(tidx2, vidx)] if cidx2 != cidx1 ]
    
                name = self.__make_name("Tv_o_coll", tidx1, cidx1, tidx2, cidx2, vidx)
                # constr = self.__problem.sum(vars) <= 1
                constr = ([1] * len(vars), vars, '<=', 1)
                self.__constraints.Tv_o_coll[(tidx1, cidx1, tidx2, cidx2, vidx)] = constr
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
        # Science targets must be either observed or go to the sink
        # If a maximum on the science target is set fo the target class, enforce that
        # in a separate constraint
        # It defaults to the total number of targets but can be overridden for each target class
        ignore_science_target_class_minimum = self.__get_debug_option('ignore_science_target_class_minimum', False)
        ignore_science_target_class_maximum = self.__get_debug_option('ignore_science_target_class_maximum', False)
        
        for target_class, vars in self.__variables.STC_o.items():
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
            min_targets = self.__target_classes[target_class].min_targets
            if not ignore_science_target_class_minimum and min_targets is not None:
                name = self.__make_name("STC_o_min", target_class)
                # constr = self.__problem.sum(vars) >= min_targets
                constr = ([1] * len(vars), vars, '>=', min_targets)
                self.__constraints.STC_o_min[target_class] = constr
                self.__add_constraint(name, constr)

            max_targets = self.__target_classes[target_class].max_targets
            if not ignore_science_target_class_maximum and max_targets is not None:
                name = self.__make_name("STC_o_max", target_class)
                # constr = self.__problem.sum(vars) <= max_targets
                constr = ([1] * len(vars), vars, '<=', max_targets)
                self.__constraints.STC_o_max[target_class] = constr
                self.__add_constraint(name, constr)

    def __create_calibration_target_class_constraints(self):
        
        
        ignore_calib_target_class_minimum = self.__get_debug_option('ignore_calib_target_class_minimum', False)
        ignore_calib_target_class_maximum = self.__get_debug_option('ignore_calib_target_class_maximum', False)
        
        for (target_class, vidx), vars in self.__variables.CTCv_o.items():
            sink = self.__variables.CTCv_sink[(target_class, vidx)]
            num_targets = len(vars)

            # The sum of all outgoing edges must be equal the number of calibration targets within each class
            name = self.__make_name("CTCv_o_sum", target_class, vidx)
            # constr = self.__problem.sum(vars + [ sink ]) == num_targets
            constr = ([1] * (len(vars) + 1), vars + [ sink ], '==', num_targets)
            self.__constraints.CTCv_o_sum[(target_class, vidx)] = constr
            self.__add_constraint(name, constr)

            # Every calibration target class must be observed a minimum number of times every visit
            min_targets = self.__target_classes[target_class].min_targets
            if not ignore_calib_target_class_minimum and min_targets is not None:
                name = self.__make_name("CTCv_o_min", target_class, vidx)
                # constr = self.__problem.sum(vars) >= min_targets
                constr = ([1] * len(vars), vars, '>=', min_targets)
                self.__constraints.CTCv_o_min[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)

            # Any calibration target class cannot be observed more tha a maximum number of times every visit
            max_targets = self.__target_classes[target_class].max_targets
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

        allow_more_visits = self.__get_netflow_option('allow_more_visits', False)

        for tidx, T_o in self.__variables.T_o.items():
            T_i = self.__variables.T_i[tidx]
            T_sink = self.__variables.T_sink[tidx]
            req_visits = self.__target_cache.req_visits[tidx]
            max_visits = len(self.__visits)

            if not allow_more_visits:
                # Require an exact number of visits
                name = self.__make_name("T_i_T_o_sum", tidx)
                # constr = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) == 0
                constr = ([ req_visits ] + [ -1 for _ in T_o ] + [ -1 ],
                          [ T_i ] + T_o + [ T_sink ], '==', 0)
                self.__constraints.T_i_T_o_sum[tidx] = constr
                self.__add_constraint(name, constr)
            else:
                # Allow for a larger number of visits than required
                name0 = self.__make_name("T_i_T_o_sum_0", tidx)
                name1 = self.__make_name("T_i_T_o_sum_1", tidx)

                # constr0 = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) <= 0
                constr0 = ([ req_visits ] + [ -1 for _ in T_o ] + [ -1 ],
                           [ T_i ] + T_o + [ T_sink ], '<=', 0)

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
            target_id = self.__target_cache.id[tidx]

            name = self.__make_name("Tv_i_Tv_o_sum", target_id, vidx)
            # constr = self.__problem.sum([ v for v in in_vars ] + [ -v[0] for v in out_vars ]) == 0
            constr = ([ 1 ] + [ -1 ] * len(out_vars), [ in_var ] + [ v for v, _ in out_vars ], '==', 0)
            self.__constraints.Tv_i_Tv_o_sum[(tidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_Cv_i_constraints(self):
        # Every cobra can observe at most one target per visit
        for (cidx, vidx), inflow in self.__variables.Cv_i.items():
            name = self.__make_name("Cv_i_sum", cidx, vidx)
            # constr = self.__problem.sum([ f for f in inflow ]) <= 1
            constr = ([1] * len(inflow), inflow, '<=', 1)
            self.__constraints.Cv_i_sum[(cidx, vidx)] = constr
            self.__add_constraint(name, constr)

    def __create_time_budget_constraints(self):
        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget

        ignore_time_budget = self.__get_debug_option('ignore_time_budget', False)
        if not ignore_time_budget:
            for budget_name, options in self.__time_budgets.items():
                budget_target_classes = set(options.target_classes)
                budget_variables = []
                # Collect all visits of targets that fall into to budget
                for (tidx, vidx), v in self.__variables.Tv_i.items():
                    target_class = self.__target_cache.target_class[tidx]
                    if target_class in budget_target_classes:
                        budget_variables.append(v)

                # TODO: This assumes a single, per visit exposure time which is fine for now
                #       but the logic could be extended further
                name = self.__make_name("Tv_i_sum", budget_name)
                # constr = self.__visit_exp_time * self.__problem.sum([ v for v in budget_variables ]) <= 3600 * options.budget
                constr = ([self.__visit_exp_time.value] * len(budget_variables), budget_variables, '<=', 3600 * options.budget)
                self.__constraints.Tv_i_sum[budget_name] = constr
                self.__add_constraint(name, constr)

    def __create_cobra_group_constraints(self, nvisits):
        # Make sure that there are enough targets in every cobra group for each visit

        ignore_cobra_group_minimum = self.__get_debug_option('ignore_cobra_group_minimum', False)
        ignore_cobra_group_maximum = self.__get_debug_option('ignore_cobra_group_maximum', False)

        for cg_name, options in self.__cobra_groups.items():
            need_min = not ignore_cobra_group_minimum and options.min_targets is not None
            need_max = not ignore_cobra_group_maximum and options.max_targets is not None
            if need_min or need_max:                
                for vidx in range(nvisits):
                    for gidx in range(options.ngroups):
                        variables = self.__variables.CG_i[(cg_name, vidx, gidx)]
                        if len(variables) > 0:
                            if need_min:
                                name = self.__make_name("Cv_CG_min", cg_name, vidx, gidx)
                                # constr = self.__problem.sum([ v for v in variables ]) >= options.min_targets
                                constr = ([1] * len(variables), variables, '>=', options.min_targets)
                                self.__constraints.CG_min[cg_name, vidx, gidx] = constr
                                self.__add_constraint(name, constr)

                            if need_max:
                                name = self.__make_name("Cv_CG_max", cg_name, vidx, gidx)
                                # constr = self.__problem.sum([ v for v in variables ]) <= options.max_targets
                                constr = ([1] * len(variables), variables, '<=', options.max_targets)
                                self.__constraints.CG_max[cg_name, vidx, gidx] = constr
                                self.__add_constraint(name, constr)

    def __create_unassigned_fiber_constraints(self, nvisits):
        # Make sure that enough fibers are kept unassigned, if this was requested
        # This is done by setting an upper limit on the sum of Cv_i edges

        ignore_reserved_fibers = self.__get_debug_option('ignore_reserved_fibers', False)
        num_reserved_fibers = self.__get_netflow_option('num_reserved_fibers', 0)

        if not ignore_reserved_fibers and num_reserved_fibers > 0:
            max_assigned_fibers = self.__bench.cobras.nCobras - num_reserved_fibers
            for vidx in range(nvisits):
                variables = [var for ((ci, vi), var) in self.__variables.Cv_i.items() if vi == vidx]
                variables = [item for sublist in variables for item in sublist]

                name = self.__make_name("Cv_i_max", vidx)
                # constr = self.__problem.sum(variables) <= max_assigned_fibers
                constr = ([1] * len(variables), variables, '<=', max_assigned_fibers)
                self.__constraints.Cv_i_max[vidx] = constr
                self.__add_constraint(name, constr)

    #endregion

    def solve(self):
        self.__problem.solve()
        self.__extract_assignments()

    def __extract_assignments(self):
        """
        Extract the fiber assignments from an LP solution
        
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
            self.__target_assignments[vidx][tidx] = cidx
            self.__cobra_assignments[vidx][cidx] = tidx

        if self.__variables is not None:
            # This only works when the netflow problem is fully built
            for (tidx, vidx), vars in self.__variables.Tv_o.items():
                for (f, cidx) in vars:
                    if self.__problem.get_value(f) > 0:
                        set_assigment(tidx, cidx, vidx)
        else:
            # This also works when only the LP problem is loaded back from a file
            # but it's much slower. It also requires that the problem is build with
            # full variable names.
            for k1, v1 in self.__problem._variables.items():
              if k1.startswith("Tv_Cv_"):
                    if self.__problem.get_value(v1) > 0:
                        _, _, tidx, cidx, vidx = k1.split("_")
                        set_assigment(int(tidx), int(cidx), int(vidx))

    # TODO: this function is redundant with the property of the same name
    #       rename to something else
    def get_cobra_assignments(self):
        """
        Return a list of arrays for each visit that contain an array with the the associated target
        index to each cobra.
        """

        return self.__cobra_assignments

    def get_target_assignments(self,  
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

        assignments : pd.DataFrame = None

        # There are 2394 cobras in total
        # There are 2604 fibers, 2394 assigned to cobras, 64 engineering fibers and 146 empty fibers
        # The design files contain 2458 rows, for the cobras plus the 64 engineering fibers

        # Internally, cidx is a 0-based index of the cobras, PFS cobraids are 1-based
        # Map the cobraIds to the corresponding fiberIds
        # These include all cobras that were part of the netflow problem

        # TODO: make sure this is updated if we leave out cobras from the netflow problem
        #       because they're broken or else

        # Should have the size of 2394
        cobraids = np.arange(self.__bench.cobras.nCobras, dtype=int) + 1
        fiberids = self.__fiber_map.cobraIdToFiberId(cobraids)

        # TODO: add fiber status, this should come from the bench config
        
        for vidx, visit in enumerate(self.__visits):
            # Unassigned cobras
            cidx = np.where(self.__cobra_assignments[vidx] == -1)[0]

            if include_unassigned_fibers:
                unassigned = pd.DataFrame(
                    {
                        'targetid': -1,
                        'pointing_idx': visit.pointing_idx,
                        'visit_idx': vidx,
                        'cobraid': cobraids[cidx],
                        'fiberid': fiberids[cidx],
                        'fp_x': np.nan,
                        'fp_y': np.nan,
                        'target_type': 'na',
                        'fiber_status': FiberStatus.GOOD
                    }
                )

                if assignments is None:
                    assignments = unassigned
                else:
                    assignments = pd.concat([ assignments, unassigned ])

            if include_engineering_fibers:
                engineering = pd.DataFrame(
                    {
                        'targetid': -1,
                        'pointing_idx': visit.pointing_idx,
                        'visit_idx': vidx,
                        'cobraid': -1,
                        'fiberid': np.where(self.__fiber_map.data['scienceFiberId'] == self.__fiber_map.ENGINEERING)[0] + 1,
                        'fp_x': np.nan,
                        'fp_y': np.nan,
                        'target_type': 'eng',
                        'fiber_status': FiberStatus.GOOD
                    }
                )

                if assignments is None:
                    assignments = engineering
                else:
                    assignments = pd.concat([ assignments, engineering ])

            # Assigned targets
            tidx = np.array([ k for k in self.__target_assignments[vidx] ], dtype=int)
            fpidx = self.__cache_to_fp_pos_map[visit.pointing_idx][tidx]

            targets = pd.DataFrame(
                {
                    'targetid': self.__target_cache.id[tidx],
                    'pointing_idx': visit.pointing_idx,
                    'visit_idx': vidx,
                    'cobraid': np.array([ cobraids[self.__target_assignments[vidx][ti]] for ti in tidx ]),
                    'fiberid': np.array([ fiberids[self.__target_assignments[vidx][ti]] for ti in tidx ]),
                    'fp_x':  self.__target_fp_pos[visit.pointing_idx][fpidx].real,
                    'fp_y':  self.__target_fp_pos[visit.pointing_idx][fpidx].imag,
                    'target_type': self.__target_cache.prefix[tidx],
                    'fiber_status': FiberStatus.GOOD
                }
            )

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
        assignments['target_type'] = assignments['target_type'].map(target_type_map)

        # Include all columns from the target list data frame, if requested
        # Convert integer columns to nullable to avoid float conversion in join
        if include_target_columns:
            assignments = assignments.join(pd_to_nullable(self.__targets), on='targetid', how='left')
            assignments = assignments.reset_index(drop=True)

        # Make sure float columns contain NaN instead of None
        pd_null_to_nan(assignments, in_place=True)

        return assignments
    
    def get_target_assignments_masks(self):
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

        all_assignments, all_fiberids = self.get_target_assignments_masks()
        num_assignments = np.sum(np.stack(all_assignments, axis=-1), axis=-1)

        targets = self.__targets.copy()
        targets['num_visits'] = num_assignments
        targets.reset_index(names='targetid', inplace=True)

        return targets

    def __get_idcol(self, catalog):
        for c in ['objid', 'skyid', 'calid']:
            if c in catalog.data.columns:
                return c
        
        raise RuntimeError()