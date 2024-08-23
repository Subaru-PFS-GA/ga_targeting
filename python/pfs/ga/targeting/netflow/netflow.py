import os
import re
import logging
from typing import Callable
from collections import defaultdict
import numpy as np
import pandas as pd
from types import SimpleNamespace  

from ics.cobraOps.Bench import Bench
from pfs.utils.coordinates.CoordTransp import CoordinateTransform

from ..util.args import *
from .pointing import Pointing
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
        * forbidden_tairs : list of list
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
    
    In the next step, the focal plane coordinates of the targets are calculated for each pointing.
    Visits of the same pointing, hence, will assume the same focal plane positions.
    
    After the setup, the ILP problem is built in `__build_ilp_problem`. This is the second most time
    consuming step of the algorithm, after the actual optimization, since it requires creating
    millions of variables and constraints.
    
    """
    
    def __init__(self, 
                 name,
                 pointings,
                 workdir=None,
                 filename_prefix=None,
                 solver_options=None,
                 netflow_options=None,
                 debug_options=None):

        # Configurable options

        self.__name = name                          # Problem name
        self.__pointings = pointings                # List of pointings
        self.__visits = None                        # List of visits
        
        self.__workdir = workdir if workdir is not None else os.getcwd()
        self.__filename_prefix = filename_prefix if filename_prefix is not None else ''
        self.__solver_options = solver_options
        self.__netflow_options = Netflow.__camel_to_snake(netflow_options)
        self.__debug_options = Netflow.__camel_to_snake(debug_options)

        # Internally used variables

        self.__bench = Bench(layout='full')

        self.__problem_type = GurobiProblem
        self.__problem = None                       # ILP problem, already wrapped

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
        self.__target_fp_pos = None                 # Target focal plane positions for each pointing
        
        self.__target_assignments = None            # List of dicts for each visit, keyed by target index
        self.__cobra_assignments = None             # List of dicts for each visit, keyed by cobra index
        self.__missed_targets = None
        self.__fully_observed_targets = None
        self.__partially_observed_targets = None

    #region Property accessors

    def __get_name(self):
        return self.__name
    
    def __set_name(self, value):
        self.__name = value

    name = property(__get_name, __set_name)

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

    def save_solution(self, filename=None):
        """Save the LP solution to a file."""

        fn = self.__get_solution_filename(filename)
        self.__problem.write_solution(fn)

    #endregion
    #region Configuration

    @staticmethod
    def __camel_to_snake(s):
        """
        If a dictionary (of settings) is provided with camelCase keys, convert them to snake_case.
        """

        if isinstance(s, dict):
            return { Netflow.__camel_to_snake(k): v for k, v in s.items() }
        elif isinstance(s, list):
            return s
        elif isinstance(s, str):
            if '_' in s:
                return s
            else:
                return re.sub(r'(?<!^)(?=[A-Z])', '_', s).lower()
        else:
            raise NotImplementedError()
        
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
        if forbidden_targets is not None:
            # Reindex targets for faster search
            df = self.__targets.set_index('id')

            for i, target_id in enumerate(forbidden_targets):            
                tidx = df.index.get_loc(target_id)
                fpp.append(tidx)
                
        return fpp
    
    def __get_forbidden_pairs_config(self):
        """
        Look up the target indices based on the target ids of forbidden pairs.
        """

        forbidden_pairs = self.__get_netflow_option('forbidden_pairs', None)
        fpp = []
        if forbidden_pairs is not None:
            # Reindex targets for faster search
            df = self.__targets.set_index('id')

            for i, pair in enumerate(forbidden_pairs):
                if len(pair) != 2:
                    raise ValueError(f"Found an incorrect number of target ids in forbidden pair list at index {i}.")
            
                tidx_list = [ df.index.get_loc(p) for p in pair ]
                fpp.append(tidx_list)
                
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

    def __calculate_fp_pos(self, pointing, ra, dec, pmra=None, pmdec=None, parallax=None, epoch=2015.5):
        cent = np.array([[ pointing.ra ], [ pointing.dec ]])
        pa = pointing.posang

        # Input coordinates. Namely. (Ra, Dec) in unit of degree for sky with shape (2, N)
        coords = np.stack([ra, dec], axis=-1)
        xyin = coords.T

        # The proper motion of the targets used for sky_pfi transformation.
        # The unit is mas/yr, the shape is (2, N)
        if pmra is not None and pmdec is not None:
            pm = np.stack([ pmra, pmdec ], axis=0)
        else:
            pm = None

        # The parallax of the coordinates used for sky_pfi transformation.
        # The unit is mas, the shape is (1, N)
        # if parallax is not None:
        #     parallax = parallax[None, :]
        # else:
        #     parallax = None

        # Observation time UTC in format of %Y-%m-%d %H:%M:%S
        obs_time = pointing.obs_time.to_value('iso')

        fp_pos = CoordinateTransform(xyin=xyin,
                                     mode="sky_pfi",
                                     # za=0.0, inr=0.0,     # These are overriden by function
                                     cent=cent,
                                     pa=pa,
                                     pm=pm,
                                     par=parallax,
                                     time=obs_time,
                                     epoch=epoch)
                
        xy = fp_pos[:2, :].T
        return xy[..., 0] + 1j * xy[..., 1]

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

        # Forbidden targets and target pairs
        self.__forbidden_targets = self.__get_forbidden_targets_config()
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
        self.__create_visits()
        self.__calculate_target_fp_pos()
        self.__build_ilp_problem()
                    
    def __calculate_exp_time(self):
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
        self.__targets['req_visits'][sci] = np.int64(np.ceil(self.__targets['exp_time'][sci] / self.__visit_exp_time))

        max_req_visits = self.__targets['req_visits'].max()
        for p in self.__pointings:
            if p.nvisits < max_req_visits:
                raise RuntimeError('Some science targets require more visits than provided.')

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
        self.__target_fp_pos = []
        for p in self.__pointings:
            fp_pos = self.__calculate_fp_pos(p, 
                                             self.__target_cache.ra,
                                             self.__target_cache.dec,
                                             # pmra=None, pmdec=None, parallax=None,
            )
            self.__target_fp_pos.append(fp_pos)

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

        logging.info("Calculating target visibilities")

        vis_targets_elbows, vis_cobras_targets = self.__calculate_visibilities()

        logging.info("Creating network topology")

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
            logging.info(f'Processing visit {visit_idx + 1}.')

            # Create calibration target class variables and define cost
            # for each calibration target class. These are different for each visit
            # because we don't necessarily need the same calibration targets at each visit.
            # > CTCv_sink_{target_class}_{visit_idx}
            logging.debug("Creating calibration target class variables.")
            self.__create_calib_target_class_variables(visit)

            # > T_Tv_{visit_idx}[target_idx]
            # > CTCv_Tv_{visit_idx}[target_idx]
            # > Tv_Cv_{visit_idx}_{cobra_idx}[target_idx]
            logging.debug("Creating target and cobra visit variables.")
            self.__create_visit_variables(visit, vis_targets_elbows[visit.pointing_idx], vis_cobras_targets[visit.pointing_idx])

            # > Tv_o_coll_{?}_{?}_{visit_idx} (endpoint collisions)
            # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx} (elbow collisions)
            logging.debug("Creating cobra collision constraints.")
            self.__create_cobra_collision_constraints(visit, vis_targets_elbows[visit.pointing_idx])

            logging.debug("Adding cobra non-allocation cost terms.")
            self.__add_cobra_non_allocation_cost(visit)

            # Add constraints for forbidden targets and forbidden pairs
            if self.__forbidden_targets is not None and len(self.__forbidden_targets) > 0:
                # > Tv_o_forb_{tidx}_{tidx}_{visit_idx}
                logging.debug("Adding forbidden target constraints")
                self.__create_forbidden_target_constraints(visit)

            if self.__forbidden_pairs is not None and len(self.__forbidden_pairs) > 0:
                # > Tv_o_forb_{tidx1}_{tidx2}_{visit_idx}
                logging.debug("Adding forbidden pair constraints")
                self.__create_forbidden_pair_constraints(visit)
        
        logging.info("Adding constraints")

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
        vis_targets_elbows = []
        vis_cobras_targets = []
        for pointing_idx, pointing in enumerate(self.__pointings):
            logging.info(f'Processing pointing {pointing_idx + 1}.')
            vte, vct = self.__get_visibility_and_elbow(self.__target_fp_pos[pointing_idx])
            vis_targets_elbows.append(vte)
            vis_cobras_targets.append(vct)

        return vis_targets_elbows, vis_cobras_targets

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
        max_visits = len(self.__visits)

        # Allow more visits than required
        T_sink =  self.__add_variable_array('T_sink', tidx, 0, max_visits)

        # Do not allow more visits than required
        # T_sink =  self.__add_variable_array('T_sink', tidx, 0, None)

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
    
    def __create_visit_variables(self, visit, vis_targets_elbows, vis_cobras_targets):
        """
        Create the variables for the visible variables for a visit.
        """

        tidx = np.array(list(vis_targets_elbows.keys()), dtype=int)

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
        for cidx, tidx_elbow in vis_cobras_targets.items():
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
            self.__variables.Tv_i[(ti, visit.visit_idx)].append(f)
            
    def __create_CTCv_Tv(self, visit, tidx):
        cost = self.__target_cache.penalty[tidx]
        name = self.__make_name('CTCv_Tv', visit.visit_idx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]
            self.__variables.Tv_i[(ti, visit.visit_idx)].append(f)
            self.__variables.CTCv_o[(target_class, visit.visit_idx)].append(f)

    def __create_Tv_Cv_CG(self, visit, cidx, tidx):
        # Calculate the cost for each target - cobra assignment
        cost = np.zeros_like(tidx, dtype=float)
        
        # Cost of a single visit
        if visit.visit_cost is not None:
            cost += visit.visit_cost

        # Cost of moving the cobra away from the center
        cobra_move_cost = self.__get_netflow_option('cobra_move_cost', None)
        if cobra_move_cost is not None:
            dist = np.abs(self.__bench.cobras.centers[cidx] - self.__target_fp_pos[visit.pointing_idx][tidx])
            cost += cobra_move_cost(dist)
        
        # Cost of closest black dots for each cobra
        if self.__black_dots is not None:
            dist = np.min(np.abs(self.__black_dots.black_dot_list[cidx] - self.__target_fp_pos[visit.pointing_idx][tidx]))
            cost += self.__black_dots.black_dot_penalty(dist)

        # Create LP variables
        name = self.__make_name("Tv_Cv", visit.visit_idx, cidx)
        vars = self.__add_variable_array(name, tidx, 0, 1, cost=cost)

        for ti in tidx:
            f = vars[ti]
            target_class = self.__target_cache.target_class[ti]

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

    def __create_cobra_collision_constraints(self, visit, vis_elbow):
        # Add constraints 
        logging.debug(f"Adding constraints for visit {visit.visit_idx}")

        # Avoid endpoint or elbow collisions
        collision_distance = self.__get_netflow_option('collision_distance', 0.0)
        elbow_collisions = self.__get_netflow_option('elbow_collisions', False)
        
        if collision_distance > 0.0:
            if not elbow_collisions:
                # > Tv_o_coll_{?}_{?}_{visit_idx}
                logging.debug("Adding endpoint collision constraints")
                self.__create_endpoint_collision_constraints(visit, vis_elbow, collision_distance)
            else:
                # > Tv_o_coll_{target_idx1}_{cobra_idx1}_{target_idx2}_{cobra_idx2}{visit_idx}
                logging.debug("Adding elbow collision constraints")
                self.__create_elbow_collision_constraints(visit, vis_elbow, collision_distance)
                
    def __create_endpoint_collision_constraints(self, visit, vis_elbow, collision_distance):
        ignore_endpoint_collisions = self.__get_debug_option('ignore_endpoint_collisions', False)

        if not ignore_endpoint_collisions:
            # Determine target indices visible by this cobra and its neighbors
            colliding_pairs = self.__get_colliding_pairs(self.__target_fp_pos[visit.pointing_idx], vis_elbow, collision_distance)
            keys = self.__variables.Tv_o.keys()
            keys = set(key[0] for key in keys if key[1] == visit.visit_idx)
            for p in colliding_pairs:
                if p[0] in keys and p[1] in keys:
                    vars = [ v for v, cidx in self.__variables.Tv_o[(p[0], visit.visit_idx)] ] + \
                           [ v for v, cidx in self.__variables.Tv_o[(p[1], visit.visit_idx)] ]
                    name = self.__make_name("Tv_o_coll", p[0], p[1], visit.visit_idx)
                    constr = self.__problem.sum(vars) <= 1
                    self.__constraints.Tv_o_coll[(p[0], p[1], visit.visit_idx)] = constr
                    self.__add_constraint(name, constr)

    def __create_elbow_collision_constraints(self, visit, vis_elbow, collision_distance):
        # Determine targets accessible by two different cobras that would cause an elbow collision
        # colliding_elbows contains a list of tidx keyed by cidx and tidx

        ignore_elbow_collisions = self.__get_debug_option('ignore_elbow_collisions', False)

        ### SLOW ###

        colliding_elbows = self.__get_colliding_elbows(self.__target_fp_pos[visit.pointing_idx], vis_elbow, collision_distance)
        for (cidx1, tidx1), tidx2_list in colliding_elbows.items():
            for f, cidx2 in self.__variables.Tv_o[(tidx1, visit.visit_idx)]:
                if cidx2 == cidx1:
                    var0 = f

            for tidx2 in tidx2_list:
                if True:  # idx2 != tidx1:
                    vars = [ var0 ]
                    vars += [ f for f, cidx2 in self.__variables.Tv_o[(tidx2, visit.visit_idx)] if cidx2 != cidx1 ]
        
                    name = self.__make_name("Tv_o_coll", tidx1, cidx1, tidx2, cidx2, visit.visit_idx)
                    constr = self.__problem.sum(vars) <= 1
                    self.__constraints.Tv_o_coll[(tidx1, cidx1, tidx2, cidx2, visit.visit_idx)] = constr
                    if not ignore_elbow_collisions:
                        self.__add_constraint(name, constr)

    def __create_forbidden_target_constraints(self, visit):
        """
        Create the constraints prohibiting individual targets being observed in any visit.
        """

        ignore_forbidden_targets = self.__get_debug_option('ignore_forbidden_targets', False)

        # All edges of visible targets, relevant for this visit
        tidx_set = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == visit.visit_idx)
        for tidx in self.__forbidden_targets:
            if tidx in tidx_set:
                flows = [ v for v, cidx in self.__variables.Tv_o[(tidx, visit.visit_idx)] ]
                name = self.__make_name("Tv_o_forb", tidx, tidx, visit.visit_idx)
                constr = self.__problem.sum(flows) == 0
                self.__constraints.Tv_o_forb[(tidx, tidx, visit.visit_idx)] = constr
                if not ignore_forbidden_targets:
                    self.__add_constraint(name, constr)

    def __create_forbidden_pair_constraints(self, visit):
        """
        Create the constraints prohibiting two targets (or individual targets) being observed
        in the same visit.
         
        This constraint is independent of the cobra collision constraints.
        """

        ignore_forbidden_pairs = self.__get_debug_option('ignore_forbidden_pairs', False)

        # All edges of visible targets, relevant for this visit
        tidx_set = set(ti for ti, vi in self.__variables.Tv_o.keys() if vi == visit.visit_idx)
        for p in self.__forbidden_pairs:
            [tidx1, tidx2] = p
            if tidx1 in tidx_set and tidx2 in tidx_set:
                flows = [ v for v, cidx in self.__variables.Tv_o[(tidx1, visit.visit_idx)] ] + \
                        [ v for v, cidx in self.__variables.Tv_o[(tidx2, visit.visit_idx)] ]
                name = self.__make_name("Tv_o_forb", tidx1, tidx2, visit.visit_idx)
                constr = self.__problem.sum(flows) <= 1
                self.__constraints.Tv_o_forb[(tidx1, tidx2, visit.visit_idx)] = constr
                if not ignore_forbidden_pairs:
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
        
        
        ignore_calib_target_class_minimum = self.__get_debug_option('ignore_calib_target_class_minimum', False)
        ignore_calib_target_class_maximum = self.__get_debug_option('ignore_calib_target_class_maximum', False)
        
        for (target_class, vidx), vars in self.__variables.CTCv_o.items():
            sink = self.__variables.CTCv_sink[(target_class, vidx)]
            num_targets = len(vars)

            # The sum of all outgoing edges must be equal the number of calibration targets within each class
            name = self.__make_name("CTCv_o_sum", target_class, vidx)
            constr = self.__problem.sum(vars + [ sink ]) == num_targets
            self.__constraints.CTCv_o_sum[(target_class, vidx)] = constr
            self.__add_constraint(name, constr)

            # Every calibration target class must be observed a minimum number of times every visit
            min_targets = self.__target_classes[target_class].min_targets
            if not ignore_calib_target_class_minimum and min_targets is not None:
                name = self.__make_name("CTCv_o_min", target_class, vidx)
                constr = self.__problem.sum(vars) >= min_targets
                self.__constraints.CTCv_o_min[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)

            # Any calibration target class cannot be observed more tha a maximum number of times every visit
            max_targets = self.__target_classes[target_class].max_targets
            if not ignore_calib_target_class_maximum and max_targets is not None:
                name = self.__make_name("CTCv_o_max", target_class, vidx)
                constr = self.__problem.sum(vars) <= max_targets
                self.__constraints.CTCv_o_max[(target_class, vidx)] = constr
                self.__add_constraint(name, constr)

    def __create_science_target_constraints(self):
        for tidx, T_o in self.__variables.T_o.items():
            T_i = self.__variables.T_i[tidx]
            T_sink = self.__variables.T_sink[tidx]
            req_visits = self.__target_cache.req_visits[tidx]
            max_visits = len(self.__visits)

            # Require an exact number of visits
            name = self.__make_name("T_i_T_o_sum", tidx)

            # Require an exact number of visits
            # constr = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] - [ T_sink ]) == 0

            # Allow for a larger number of visits than required
            constr0 = self.__problem.sum([ req_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) <= 0
            constr1 = self.__problem.sum([ max_visits * T_i ] + [ -v for v in T_o ] + [ -T_sink ]) >= 0

            self.__constraints.T_i_T_o_sum[(tidx, 0)] = constr0
            self.__constraints.T_i_T_o_sum[(tidx, 1)] = constr1
            self.__add_constraint(name, constr0)
            self.__add_constraint(name, constr1)

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
                constr = self.__visit_exp_time * self.__problem.sum([ v for v in budget_variables ]) <= 3600 * options.budget
                self.__constraints.Tv_i_sum[budget_name] = constr
                self.__add_constraint(name, constr)

    def __create_cobra_group_constraints(self, nvisits):
        # Make sure that there are enough targets in every cobra group for each visit

        ignore_cobra_group_minimum = self.__get_debug_option('ignore_cobra_group_minimum', False)
        ignore_cobra_group_maximum = self.__get_debug_option('ignore_cobra_group_maximum', False)

        for name, options in self.__cobra_groups.items():
            need_min = not ignore_cobra_group_minimum and options.min_targets is not None
            need_max = not ignore_cobra_group_maximum and options.max_targets is not None
            if need_min or need_max:                
                for vidx in range(nvisits):
                    for gidx in range(options.ngroups):
                        variables = self.__variables.CG_i[(name, vidx, gidx)]
                        if len(variables) > 0:
                            if need_min:
                                name = self.__make_name("Cv_CG_min", name, vidx, gidx)
                                constr = self.__problem.sum([ v for v in variables ]) >= options.min_targets
                                self.__constraints.CG_min[name, vidx, gidx] = constr
                                self.__add_constraint(name, constr)

                            if need_max:
                                name = self.__make_name("Cv_CG_max", name, vidx, gidx)
                                constr = self.__problem.sum([ v for v in variables ]) <= options.max_targets
                                self.__constraints.CG_max[name, vidx, gidx] = constr
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
                constr = self.__problem.sum(variables) <= max_assigned_fibers
                self.__constraints.Cv_i_max[vidx] = constr
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
        """
        Find science targets for each science target class that have been missed.
        These can be identified by the T_i variables that are zero.
        """

        self.__missed_targets = { k: [] for k, options in self.__target_classes.items() if options.prefix == 'sci' }

        if self.__variables is not None:
            for tidx, f in self.__variables.T_i.items():
                if self.__problem.get_value(f) == 0:
                    target_class = self.__target_cache.target_class[tidx]
                    self.__missed_targets[target_class].append(tidx)
        else:
            raise NotImplementedError()
        
    def __extract_fully_observed_science_targets(self):
        """
        Find science targets that are allocated to fibers during at least the number of visits
        as requested.
        """

        self.__fully_observed_targets = { k: [] for k, options in self.__target_classes.items() if options.prefix == 'sci' }

        if self.__variables is not None:
            for tidx, f in self.__variables.T_sink.items():
                g = self.__variables.T_i[tidx]
                if self.__problem.get_value(f) == 0 and self.__problem.get_value(g) == 1:
                    target_class = self.__target_cache.target_class[tidx]
                    self.__fully_observed_targets[target_class].append(tidx)
        else:
            raise NotImplementedError()
        
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
    
    # TODO: this function is redundant with the property of the same name
    #       rename to something else
    def get_cobra_assignments(self):
        """
        Return a list of arrays for each visit that contain an array with the the associated target
        index to each cobra.
        """

        return self.__cobra_assignments
    
    # TODO: this function is redundant with the property of the same name
    #       rename to something else
    def get_target_assignments(self):
        """
        Returns a list of data frames, for each visit, of the target assignments including
        positions, fiberid and other information.
        """

        assignments = []
        for i, visit in enumerate(self.__visits):
            # Target indices
            idx = np.array([ k for k in self.__target_assignments[i] ], dtype=int)

            targets = pd.DataFrame(
                    {
                        f'targetid': [ np.int64(self.__targets.iloc[ti]['id']) for ti in idx ],
                        'fiberid': np.array([ self.__target_assignments[i][k] for k in idx ]),
                        'fp_x':  self.__target_fp_pos[visit.pointing_idx][idx].real,
                        'fp_y':  self.__target_fp_pos[visit.pointing_idx][idx].imag,
                    }
                ).astype({ 'fiberid': pd.Int32Dtype() }).set_index(f'targetid')

            # Find the index of each object in the catalog and return the matching fiberid
            assign = self.__targets.set_index('id').join(targets, how='inner')
            assign = assign.reset_index(names=['targetid'])
            assignments.append(assign)

        return assignments
    
    def get_target_assignments_masks(self):
        """
        Returns a mask that indexes the assigned targets within the target list.
        """

        masks = []
        fiberids = []
        for i, p in enumerate(self.__visits):
            m = np.full(len(self.__targets), False, dtype=bool)
            m[list(self.__target_assignments[i].keys())] = True
            masks.append(m)
            fiberids.append(np.array(self.__target_assignments[i].values()))

        return masks, fiberids
    
    def get_catalog_assignments(self, catalog):
        """
        Cross-references the provided catalog with the target assignments and returns a list of
        data frames for each visit with two columns: the object IDs and the fiber ID.

        Parameters:
        -----------
        catalog: Catalog
            The catalog to cross-reference with the target assignments.
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

    def get_catalog_assignments_masks(self, catalog):
        """
        Returns a mask that indexes the selected targets within a catalog.
        """

        idcol = self.__get_idcol(catalog)
        assignments = self.get_catalog_assignments(catalog)
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