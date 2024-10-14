import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from test_base import TestBase

from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.Bench import Bench
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID

import pfs.ga.targeting
from pfs.ga.targeting.instrument import SubaruPFI
from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.netflow import Netflow
from pfs.ga.targeting.data import Observation

class NetflowTest(TestBase):
    def test_check_pointing_visibility(self):
        # Ursa Minor dSph visible
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        nf = Netflow(f'test', instrument, [ pointing ])
        nf._Netflow__check_pointing_visibility()

        # Below horizon
        pointing = Pointing(0, 0, posang=0, obs_time=Time("2016-04-03T08:00:00Z"))
        nf = Netflow(f'test', instrument, [ pointing ])
        with self.assertRaises(AssertionError):
            nf._Netflow__check_pointing_visibility()

    def test_filter_targets(self):
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1)
        nf = Netflow(f'test', instrument, [ pointing ])

        obs = self.load_test_observation()
        nf._Netflow__filter_targets(obs.data, [ pointing])

    def test_cache_targets(self):
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1, exp_time=1200)
        obs = self.load_test_observation()
        nf = Netflow(f'test', instrument, [ pointing ])
        nf.append_science_targets(obs, exp_time=1200, priority=1)
        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        
        nf._Netflow__cache_targets()

    def test_calculate_target_fp_pos(self):
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1, exp_time=1200)
        obs = self.load_test_observation()
        nf = Netflow(f'test', instrument, [ pointing ])
        nf.append_science_targets(obs, exp_time=1200, priority=1)
        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        nf._Netflow__cache_targets()

        nf._Netflow__calculate_target_fp_pos()
        nf._Netflow__check_target_fp_pos()

    def get_gurobi_options(self):
        return dict(
            seed=0,                 # random seed
            presolve=2,             # agressiveness of presolve which tries to eliminate variables from the LP problem
            method=3,               # 3 means concurrent, 4 means deterministic concurrent
            degenmoves=0,           # degenerate simplex moves, set to 0 to prevent too much time to be spent on trying to improve the current solution
            heuristics=0.5,         # how much of the time to spend by performing heuristics
            mipfocus=1,             # mipfocus=1 is balanced toward finding more feasible solutions
                                    # mipfocus=2 is balanced toward proving that the current solution is the best
                                    # mipfocus=3 is to be used when the objection bound is moving very slowly
            mipgap=0.01,            # relative stopping criterion for bounds on the objective
            LogToConsole=1,         # 
            timelimit=300           # in sec
        )
    
    def get_target_classes(self):
        target_classes = {
            # 'sky': dict(
            #     prefix = 'sky',
            #     min_targets = 240,
            #     max_targets = 320,
            #     non_observation_cost = 0,
            # ),
            # 'cal': dict(
            #     prefix = 'cal',
            #     min_targets = 40,
            #     max_targets = 240,
            #     non_observation_cost = 0,
            #     calib = True,
            # ),
        }
    
        for i in range(10):
            target_classes[f'sci_P{i}'] = dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = max(20 - 2 * i, 1),
                partial_observation_cost = 1e5,
            )

        target_classes[f'sci_P0']['non_observation_cost'] = 1000
        target_classes[f'sci_P1']['non_observation_cost'] = 500
        target_classes[f'sci_P2']['non_observation_cost'] = 200
        target_classes[f'sci_P3']['non_observation_cost'] = 100
        target_classes[f'sci_P4']['non_observation_cost'] = 100
        target_classes[f'sci_P5']['non_observation_cost'] = 100
        target_classes[f'sci_P6']['non_observation_cost'] = 100
        target_classes[f'sci_P7']['non_observation_cost'] = 50
        target_classes[f'sci_P8']['non_observation_cost'] = 10
        target_classes[f'sci_P9']['non_observation_cost'] = 0

        return target_classes

    def get_cobra_groups(self):
        cobra_groups = {
            # 'location': dict(
            #     groups = np.random.randint(4, size=ncobras),
            #     target_classes = [ 'sky' ],
            #     min_targets = 40,
            #     max_targets = 80,
            #     non_observation_cost = 10,
            # ),
            # 'instrument': dict(
            #     groups = np.random.randint(4, size=ncobras),
            #     target_classes = [ 'sky' ],
            #     min_targets = 10,
            #     max_targets = 25,
            #     non_observation_cost = 10,
            # )
        }

        return cobra_groups
    
    def get_netflow_options(self, target_classes, cobra_groups):
        netflow_options = dict(
            # Add a penalty if the target is too close to a black dot
            black_dot_penalty = lambda dist: 1.0,

            # Penalize cobra moves with respect to cobra center
            cobra_move_cost = lambda dist: 1.0,

            fiber_non_allocation_cost = 1e5,

            collision_distance = 2.0,
            elbow_collisions = True,
            # forbidden_targets = [
            #     43218108431
            # ],
            # forbidden_pairs = [
            #     [43486543901, 43218108431],
            # ],

            target_classes = target_classes,
            cobra_groups = cobra_groups,

            time_budgets = {
                'science': dict(
                    target_classes = [ 'sci_P0', 'sci_P1', 'sci_p2' ],
                    budget = 5  # hr
                )
            },

            num_reserved_fibers = 0,

            # This will only be used when netflow is rewritten to step-by-step fiber assignment
            # constrain_already_observed = False,

            # Allow more visits than minimally required
            allow_more_visits = True,

            epoch = 2016, # all catalogs must match
            ignore_proper_motion = True,

            # FPI configuration
            fiberids_path = os.path.join(os.path.dirname(pfs.utils.__file__), '../../../data/fiberids')
        )

        return netflow_options
    
    def get_debug_options(self):
        debug_options = dict(
            ignoreEndpointCollisions = False,
            ignoreElbowCollisions = False,
            ignoreForbiddenPairs = False,
            ignoreForbiddenSingles = False,
            ignoreCalibTargetClassMinimum = False,
            ignoreCalibTargetClassMaximum = False,
            ignoreScienceTargetClassMinimum = False,
            ignoreScienceTargetClassMaximum = False,
            ignoreTimeBudget = False,
            ignoreCobraGroupMinimum = False,
            ignoreCobraGroupMaximum = False,
            ignoreReservedFibers = False,
        )

        return debug_options

    def test_solve(self):
        # This is a full coverage test, does not validate the results and
        # correctness of the individual features.

        instrument_options = {}
        gurobi_options = self.get_gurobi_options()
        target_classes = self.get_target_classes()
        cobra_groups = self.get_cobra_groups()
        netflow_options = self.get_netflow_options(target_classes, cobra_groups)
        debug_options = self.get_debug_options()

        instrument = SubaruPFI(instrument_options=instrument_options)
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=2, exp_time=1200)
        nf = Netflow(f'test', instrument, [ pointing ],
                     netflow_options=netflow_options,
                     solver_options=gurobi_options,
                     debug_options=debug_options)

        obs = self.load_test_observation()
        nf.append_science_targets(obs, exp_time=1200, priority=1)

        # Ignore some targets, this is to test the target index mappings
        # This is for testing, do not access the target cache directly

        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        nf._Netflow__cache_targets()

        # Pick the first and last target that's inside the pointings
        id1, id2 = nf._Netflow__target_cache.id[0], nf._Netflow__target_cache.id[-1]
        print(id1, id2)

        # Pick the first and last target
        id3, id4 = nf._Netflow__targets.index[0], nf._Netflow__targets.index[-1]
        print(id3, id4)

        nf.netflow_options['forbidden_targets'] = [ id1, id2, id3, id4 ]
        nf.netflow_options['forbidden_pairs'] = [ [ id1, id2 ], [ id3, id4 ] ]

        nf._Netflow__targets.reset_index(names='id', inplace=True)
        
        nf.build()
        nf.solve()

        assignments = nf.get_target_assignments(include_target_columns=True,
                                                include_unassigned_fibers=True,
                                                include_engineering_fibers=True)
        
        summary = nf.get_target_assignment_summary()