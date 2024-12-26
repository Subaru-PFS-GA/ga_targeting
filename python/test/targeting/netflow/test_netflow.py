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
from pfs.ga.targeting.config import NetflowConfig

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

        # Bootes I is visible at a very large airmass
        pointing = Pointing(210.025, 14.5, posang=30, obs_time=Time("2025-01-24T11:00:00Z"))
        nf = Netflow(f'test', instrument, [ pointing ])
        with self.assertRaises(AssertionError):
            nf._Netflow__check_pointing_visibility()

    def test_filter_targets(self):
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 14.5, posang=30, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1)
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
        pointing = Pointing(227.1, 67.25, posang=30, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1, exp_time=1200)
        obs = self.load_test_observation()
        config = NetflowConfig.default()
        nf = Netflow(f'test', instrument, [ pointing ],
                     netflow_options=config.netflow_options,
                     debug_options=config.debug_options)
        nf.append_science_targets(obs, exp_time=1200, priority=1)
        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        nf._Netflow__cache_targets()

        nf._Netflow__calculate_target_fp_pos()
        nf._Netflow__check_target_fp_pos()

    def test_calculate_target_fp_pos_error(self):
        instrument = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1, exp_time=1200)
        obs = self.load_test_observation()
        config = NetflowConfig.default()
        nf = Netflow(f'test', instrument, [ pointing ],
                     netflow_options=config.netflow_options,
                     debug_options=config.debug_options)
        nf.append_science_targets(obs, exp_time=1200, priority=1)
        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        nf._Netflow__cache_targets()

        nf._Netflow__calculate_target_fp_pos()

        # Focal plane position is set such that observations don't cover it completely
        self.assertRaises(AssertionError, nf._Netflow__check_target_fp_pos)
        
    def test_solve(self):
        # This is a full coverage test, does not validate the results and
        # correctness of the individual features.

        config = NetflowConfig.default()
        config.load(os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../configs/netflow/example.py'))

        instrument = SubaruPFI(instrument_options=config.instrument_options)
        pointing = Pointing(227.1, 67.25, posang=30, obs_time=Time("2024-06-10T00:00:00.0Z"), nvisits=1, exp_time=1200)
        nf = Netflow(f'test', instrument, [ pointing ],
                     netflow_options=config.netflow_options,
                     solver_options=config.gurobi_options,
                     debug_options=config.debug_options)

        obs = self.load_test_observation()
        nf.append_science_targets(obs, exp_time=1200)

        # Ignore some targets, this is to test the target index mappings
        # This is for testing, do not access the target cache directly

        nf._Netflow__calculate_exp_time()
        nf._Netflow__calculate_target_visits()
        nf._Netflow__cache_targets()

        # Pick the first and last target that's inside the pointings
        (id1, id2, id3, id4) = nf._Netflow__target_cache.id[[0, 1, -2, -1]]

        nf.netflow_options.forbidden_targets = [ id1, id2, id3, id4 ]
        nf.netflow_options.forbidden_pairs = [ [ id1, id2 ], [ id3, id4 ] ]

        nf._Netflow__targets.reset_index(names='targetid', inplace=True)
        
        nf.build()
        nf.solve()

        assignments = nf.get_target_assignments(include_target_columns=True,
                                                include_unassigned_fibers=True,
                                                include_engineering_fibers=True)
        
        summary = nf.get_target_assignment_summary()