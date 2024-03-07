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
                 netflow_options=None):

        # Configurable options

        self.__name = name                          # Problem name
        self.__visits = visits                      # List of pointings
        self.__target_classes = {}
        self.__target_forbidden_pairs = None
        
        self.__workdir = workdir if workdir is not None else os.getcwd()
        self.__filename_prefix = filename_prefix if filename_prefix is not None else ''
        self.__solver_options = solver_options
        self.__netflow_options = netflow_options

        # Internally used variables

        self.__telescope_type = SubaruWFC
        self.__positioner_type = SubaruPFI
        self.__positioner = None

        self.__problem_type = GurobiProblem
        self.__problem = None                       # ILP problem, already wrapped

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
        if key in self.__netflow_options:
            return self.__netflow_options[key]
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
        
    def __append_targets(self, catalog, id_column, prefix, exp_time=None, priority=None, penalty=None, mask=None, filter=None):
        df = catalog.get_data(mask=mask, filter=filter)

        targets = pd.DataFrame({ 'id': df[id_column],
                                 'RA': df['RA'],
                                 'Dec': df['Dec'] })
        
        # TODO: add pm, epoch

        if exp_time is not None:
            targets['exp_time'] = exp_time
        else:
            targets['exp_time'] = df['exp_time']

        if priority is not None:
            targets['priority'] = priority
        else:
            targets['priority'] = df['priority']

        # This is a per-object penalty for observing calibration targets
        if penalty is not None:
            targets['penalty'] = penalty
        elif 'penalty' in df.columns:
            targets['penalty'] = df['penalty']
        else:
            targets['penalty'] = 0

        targets['prefix'] = prefix

        if prefix == 'sci':
            targets['class'] = targets[['prefix', 'priority']].apply(lambda r: f"{r['prefix']}_P{r['priority']}", axis=1)
        else:
            targets['class'] = prefix

        if self.__targets is None:
            self.__targets = targets
        else:
            self.__targets = pd.concat(self.__targets, targets)
            
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

        self.__positioner = self.__positioner_type()

        # Validation steps
        self.__validate_target_classes()

        # Build the problem
        self.__calculate_target_fp_pos()
        self.__collect_target_forbidden_pairs()
        self.__calculate_exp_time()
        self.__calculate_target_visits()

        self.__build_ilp_problem()
        
    def __validate_target_classes(self):
        # sanity check for science targets: make sure that partialObservationCost
        # is larger than nonObservationCost
        for key, val in self.__target_classes.items():
            if not val["calib"] and val["partialObservationCost"] < val["nonObservationCost"]:
                raise ValueError(
                    "Found a target class where partialObservationCost "
                    "is smaller than nonObservationCost")

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

            ra = np.array(self.__targets['RA'])
            dec = np.array(self.__targets['Dec'])

            # Calculate focal plane positions of targets
            fpp, _ = telescope.world_to_fp_pos(ra, dec)
            self.__target_fp_pos.append(fpp[..., 0] + 1j * fpp[..., 1])

    def __collect_target_forbidden_pairs(self):
        self.__target_forbidden_pairs = []
        for i, pointing in enumerate(self.__visits):
            self.__target_forbidden_pairs.append([])

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

    def __make_name(self, *parts):
        return "_".join([str(x) for x in parts])

    def __build_ilp_problem(self):

        # Cobra movement cost
        # TODO: Maybe use this to prevent cobra from moving?
        def cobra_move_cost(dist):
            return 0.0 * dist   # Do not penalize cobra moves.

        # Number of visits
        nvisits = len(self.__visits)

        self.__problem = self.__problem_type(self.__name, self.__solver_options)

        Cv_i = defaultdict(list)        # Cobra visit inflows, key: (cobra_idx, visit_idx)
        Tv_o = defaultdict(list)        # Target visit outflows
        Tv_i = defaultdict(list)        # Target visit inflows, key: (target_idx, visit_idx)
        T_o = defaultdict(list)         # Target outflows (only science targets)
        T_i = defaultdict(list)         # Target inflows (only science targets), key: (target_idx)
        CTCv_o = defaultdict(list)      # Calibration target class visit outflows, key: (target_class, visit_idx)
        STC_o = defaultdict(list)       # Science Target outflows, key: (target_class)

        logging.debug("Creating network topology")

        # Generate sinks for cobra location groups used to distribute sky fibers
        # TODO: since these work similarly, we can just merge them as cobra_groups and have any number of them
        cobra_location_group = self.__get_netflow_option('cobra_location_group', None)
        cobra_location_group_max, cobra_location_group_variables = self.__add_cobra_group_sinks(nvisits,
                                                                cobra_location_group,
                                                                'locationGroupPenalty',
                                                                'locgroup_sink')

        # Generate sinks for cobra instrument groups used to distribute sky fibers
        cobra_instrument_region = self.__get_netflow_option('cobra_instrument_region', None)
        cobra_instrument_region_max, cobra_instrument_region_variables = self.__add_cobra_group_sinks(nvisits,
                                                                cobra_instrument_region,
                                                                'instrumentRegionPenalty',
                                                                'instregion_sink')
        
        # Closest black dots for each cobra
        black_dot_penalty = self.__get_netflow_option('blackDotPenalty', None)
        if black_dot_penalty is not None:
            black_dot_list = self.__positioner.nf_get_closest_dots()
        else:
            black_dot_list = None
        
        # Create nodes for every visit and calibration target class
        for key, value in self.__target_classes.items():
            if value["calib"]:
                for ivis in range(nvisits):
                    f = self.__problem.add_variable(self.__make_name("CTCv_sink", key, ivis), 0, None)
                    CTCv_o[(key, ivis)].append(f)
                    self.__problem.add_cost(f * value["nonObservationCost"])

        # Create the rest of the variables
        for ivis in range(nvisits):
            logging.debug(f"  exposure {ivis + 1}")
            logging.debug("Calculating visibilities")

            # TODO: rename
            vis_elbow = self.__positioner.nf_get_visibility_and_elbow(self.__target_fp_pos[ivis])
            for tidx, cidx_elbow in vis_elbow.items():
                target = self.__targets.iloc[tidx]
                target_id = target['id']
                target_prefix = target['prefix']
                tc = target['class']
                target_penalty = target['penalty']
                target_done_visits = target['done_visits']

                if target_prefix == 'sci':
                    # Target node to target visit node
                    f = self.__problem.add_variable(self.__make_name("T_Tv", target_id, ivis), 0, 1)
                    T_o[tidx].append(f)
                    Tv_i[(tidx, ivis)].append(f)

                    # TODO: why is it inside a loop but in an if?
                    if len(T_o[tidx]) == 1:  # freshly created

                        # Science Target class node to target node
                        # TODO: This is a bit fishy here. If the target is marked as already observed
                        #       by done_visits > 0, why is a new variable added at all?
                        lo = 1 if target_done_visits > 0 else 0
                        f = self.__problem.add_variable(self.__make_name("STC_T", tc, target_id), lo, 1)
                        T_i[tidx].append(f)
                        STC_o[tc].append(f)

                        if len(STC_o[tc]) == 1:  # freshly created
                            # Science Target class node to sink
                            f = self.__problem.add_variable(self.__make_name("STC_sink", tc), 0, None)
                            STC_o[tc].append(f)

                            self.__problem.add_cost(f * self.__target_classes[tc]["nonObservationCost"])
                        
                        # Science Target node to sink
                        f = self.__problem.add_variable(self.__make_name("ST_sink", target_id), 0, None)
                        T_o[tidx].append(f)

                        self.__problem.add_cost(f * self.__target_classes[tc]["partialObservationCost"])

                elif target_prefix in ['sky', 'cal']:
                    # Calibration Target class node to target visit node
                    f = self.__problem.add_variable(self.__make_name("CTCv_Tv", tc, target_id, ivis), 0, 1)
                    Tv_i[(tidx, ivis)].append(f)
                    CTCv_o[(tc, ivis)].append(f)
                    
                    if target_penalty != 0:
                        self.__problem.add_cost(f * target_penalty)
                else:
                    raise NotImplementedError()

                # Target visit node to cobra visit node
                for (cidx, _) in cidx_elbow:
                    f = self.__problem.add_variable(self.__make_name("Tv_Cv", tidx, cidx, ivis), 0, 1)
                    Cv_i[(cidx, ivis)].append(f)
                    Tv_o[(tidx, ivis)].append((f, cidx))

                    if cobra_location_group is not None and target_prefix == "sky":
                        cobra_location_group_variables[ivis][cobra_location_group[cidx]].append(f)

                    if cobra_instrument_region is not None and target_prefix == "sky":
                        cobra_instrument_region_variables[ivis][cobra_instrument_region[cidx]].append(f)

                    # Cost of the visit
                    tcost = self.__visits[ivis].obs_cost

                    # Cost of moving the cobra
                    if cobra_move_cost is not None:
                        dist = np.abs(self.__positioner.bench.cobras.centers[cidx] - self.__target_fp_pos[ivis][tidx])
                        tcost += cobra_move_cost(dist)
                    
                    # TODO: either zero tcost out of skip this
                    # if tcost != 0:
                    #     self.__problem.add_cost(f * tcost)

                    # Cost of closest black dots for each cobra
                    if black_dot_penalty is not None:
                        dist = np.min(np.abs(black_dot_list[cidx] - self.__target_fp_pos[ivis][tidx]))
                        tcost += black_dot_penalty(dist)
                    
                    if tcost != 0:
                        self.__problem.add_cost(f * tcost)

            # If requested, penalize non-allocated fibers
            fiber_non_allocation_cost = self.__get_netflow_option('fiberNonAllocationCost', 0)
            if fiber_non_allocation_cost != 0:
                # TODO: consider storing Cv_i organized by visit instead of cobra as well
                relevant_vars = [ var for ((ci, vi), var) in Cv_i.items() if vi == ivis ]
                relevant_vars = [ item for sublist in relevant_vars for item in sublist ]
                self.__problem.add_cost(fiber_non_allocation_cost * (self.__positioner.bench.cobras.nCobras - self.__problem.sum(relevant_vars)))

            # Add constraints 
            logging.debug(f"Adding constraints for visit {ivis}")

            collision_constraints = []

            # Avoid endpoint collisions
            collision_distance = self.__get_netflow_option('collision_distance', 0.0)
            elbow_collisions = self.__get_netflow_option('elbow_collisions', False)
            
            if collision_distance > 0.0:
                logging.debug("Adding collision constraints")
                
                if not elbow_collisions:
                    # Determine target indices visible by this cobra and its neighbors
                    colliding_pairs = self.__positioner.nf_get_colliding_pairs(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
                    keys = Tv_o.keys()
                    keys = set(key[0] for key in keys if key[1] == ivis)
                    for p in colliding_pairs:
                        if p[0] in keys and p[1] in keys:
                            flows = [v[0] for v in Tv_o[(p[0], ivis)] + Tv_o[(p[1], ivis)]]
                            tname0 = self.__targets.iloc[p[0]]['id']
                            tname1 = self.__targets.iloc[p[1]]['id']
                            collision_constraints.append([
                                self.__make_name("Coll_", tname0, tname1, ivis),
                                self.__problem.sum(flows) <= 1
                            ])
                else:
                    colliding_elbows = self.__positioner.nf_get_colliding_elbows(self.__target_fp_pos[ivis], vis_elbow, collision_distance)
                    for (cidx, tidx1), tidx2 in colliding_elbows.items():
                        for f2, cidx2 in Tv_o[(tidx1, ivis)]:
                            if cidx2 == cidx:
                                flow0 = f2
                        for idx2 in tidx2:
                            if True:  # idx2 != tidx1:
                                flows = [ flow0 ]
                                flows += [ f2 for f2, cidx2 in Tv_o[(idx2, ivis)] if cidx2 != cidx ]
                                collision_constraints.append([
                                    self.__make_name("Coll_", cidx, cidx2, ivis),
                                    self.__problem.sum(flows) <= 1 ])

            # add constraints for forbidden pairs of targets
            if self.__target_forbidden_pairs is not None:
                logging.debug("Adding forbidden pair constraints")

                keys = Tv_o.keys()
                keys = set(key[0] for key in keys if key[1] == ivis)
                for p in self.__target_forbidden_pairs[ivis]:
                    if len(p) == 2:
                        if p[0] in keys and p[1] in keys:
                            flows = [v[0] for v in Tv_o[(p[0], ivis)] + Tv_o[(p[1], ivis)]]
                            tname0 = self.__targets.iloc[p[0]]['id']
                            tname1 = self.__targets.iloc[p[1]]['id']
                            collision_constraints.append([
                                self._make_name("forbiddenPair_", tname0, tname1, ivis),
                                self.__problem.sum(flows) <= 1
                            ])
                    elif len(p) == 1:
                        if p[0] in keys:
                            flows = [v[0] for v in Tv_o[(p[0], ivis)]]
                            tname0 = self.__targets.iloc[p[0]]['id']
                            collision_constraints.append([
                                self._make_name("forbiddenPair_", tname0, ivis),
                                self.__problem.sum(flows) == 0
                            ])
                    else:
                        raise RuntimeError("oops")
                    
        for (name, constr) in collision_constraints:
            # # We add the collision constraints as lazy in the hope that this will
            # # speed up the solution
            # self.__problem.add_lazy_constraint(c[0], c[1])

            self.__problem.add_constraint(name, constr)

        # Every cobra can observe at most one target per visit
        for (cidx, vidx), inflow in Cv_i.items():
            name = self.__make_name("Cvlim", cidx, vidx)
            constr = self.__problem.sum([ f for f in inflow ]) <= 1
            self.__problem.add_constraint(name, constr)
            
        # Every calibration target class must be observed a minimum number of times
        # every visit
        for (tc, vidx), vars in CTCv_o.items():
            name = self.__make_name("CTCv_min", tc, vidx)
            constr = self.__problem.sum([ v for v in vars ]) >= self.__target_classes[tc]["numRequired"]
            self.__problem.add_constraint(name, constr)

        # Inflow and outflow at every Tv node must be balanced
        for (tidx, vidx), ivars in Tv_i.items():
            ovars = Tv_o[(tidx, vidx)]
            name = self.__make_name("TvIO", self.__targets.iloc[tidx]['id'], vidx)
            constr = self.__problem.sum([ v for v in ivars ] + [ -v[0] for v in ovars ]) == 0
            self.__problem.add_constraint(name, constr)

        # Inflow and outflow at every T node must be balanced
        for tidx, ivars in T_i.items():
            ovars = T_o[tidx]
            nvis = max(0, self.__targets.iloc[tidx]['req_visits'] - self.__targets.iloc[tidx]['done_visits'])
            name = self.__make_name("TIO", self.__targets.iloc[tidx]['id'])
            constr = self.__problem.sum([ nvis * v for v in ivars ] + [ -v for v in ovars ]) == 0
            self.__problem.add_constraint(name, constr)

        # The maximum number of targets to be observed within a science target class
        # It defaults to the total number of targets but can be overridden for each target class
        for tc, vars in STC_o.items():
            n_obs = len(vars) - 1
            if "nobs_max" in self.__target_classes[tc]:
                n_obs = self.__target_classes[tc]["nobs_max"]
            name = self.__make_name("ST", key[0], key[1])
            constr = self.__problem.sum([ v for v in vars ]) == n_obs
            self.__problem.add_constraint(name, constr)

        # Science targets inside a given program must not get more observation time
        # than allowed by the time budget
        obsprog_time_budget = self.__get_netflow_option('obsprog_time_budget', None)
        if obsprog_time_budget is not None:
            budget_variables = defaultdict(list)
            for (tidx, ivis), val in Tv_i.items():
                tc = self.__targets.iloc[tidx]['class']
                for classes in obsprog_time_budget.keys():
                    if tc in classes:
                        budget_variables[classes].append(val[0]) 

            # TODO: This assumes a single per visit exposure time which is fine for now
            #       but it could be extended further
            budget_count = 0
            for classes, val in budget_variables.items():
                name = self.__make_name("budget", budget_count)
                constr = self.__visit_exp_time * self.__problem.sum([ v for v in val ]) <= 3600 * obsprog_time_budget[classes]
                self.__problem.add_constraint(name, constr)
                budget_count += 1

        # Make sure that there are enough sky targets in every Cobra location group
        minSkyTargetsPerLocation = self.__get_netflow_option('minSkyTargetsPerLocation', 0)
        if cobra_location_group is not None and len(cobra_location_group) > 0:
            for ivis in range(nvisits):
                for i in range(cobra_location_group_max + 1):
                    name = self.__make_name("LocGrp", ivis, i)
                    constr = self.__problem.sum([ v for v in cobra_location_group_variables[ivis][i] ]) >= minSkyTargetsPerLocation
                    self.__problem.add_constraint(name, constr)

        # Make sure that there are enough sky targets in every Cobra instrument region
        minSkyTargetsPerInstrumentRegion = self.__get_netflow_option('minSkyTargetsPerInstrumentRegion', 0)
        if cobra_instrument_region is not None and len(cobra_instrument_region) > 0:
            for ivis in range(nvisits):
                for i in range(cobra_instrument_region_max + 1):
                    name = self.__make_name("InstReg", ivis, i)
                    constr = self.__problem.sum([ v for v in cobra_instrument_region_variables[ivis][i]] ) >= minSkyTargetsPerInstrumentRegion
                    self.__problem.add_constraint(name, constr)

        # Make sure that enough fibers are kept unassigned, if this was requested
        numReservedFibers = self.__get_netflow_option('numReservedFibers', 0)      
        if numReservedFibers > 0:
            maxAssignableFibers = self.__positioner.bench.cobras.nCobras - numReservedFibers
            for ivis in range(nvisits):
                relevantVars = [var for (keys, var) in Cv_i.items() if keys[1] == ivis]
                relevantVars = [item for sublist in relevantVars for item in sublist]
                name = self.__make_name("FiberLimit", ivis)
                constr = self.__problem.sum(relevantVars) <= maxAssignableFibers
                self.__problem.add_constraint(name, constr)

        pass

    def __add_cobra_group_sinks(self, nvisits, group, penalty_option_key, prefix):
        if group is not None and len(group) > 0:
            max_group = max(group)
            variables = [[[] for _ in range(max_group + 1)] for _ in range(nvisits)]
            # Add overflow arcs.
            for i in range(nvisits):
                for j in range(max_group + 1):
                    f = self.__problem.add_variable(self.__make_name(prefix, i, j), 0, None)
                    variables[i][j].append(f)
                    self.__problem.add_cost(f * self.__get_netflow_option(penalty_option_key, 20))
        else:
            max_group = None
            variables = None

        return max_group, variables

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