import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from ics.cobraOps.TargetGroup import TargetGroup

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruWFC, SubaruPFI
from pfs.ga.targeting.netflow import CollisionSimulator

from test_base import TestBase

class CollisionSimulatorTest(TestBase):
    def get_instrument(self):
        # Load the instrument calibration data

        p = Pointing(210.025, 14.5, posang=30, obs_time=datetime(2025, 1, 25, 13, 30, 0))
        wfc = SubaruWFC(p)
        pfi = SubaruPFI(wfc, instrument_options={ 'layout': 'calibration' })

        return pfi

    def generate_targets(self, pfi):
        # Generate random focal plane positions

        # For each cobra, generate a random theta and phi
        theta = np.random.rand(pfi.bench.cobras.nCobras) * 2 * np.pi
        phi = np.random.rand(pfi.bench.cobras.nCobras) * 2 * np.pi - np.pi

        # Calculate the elbow and fiber focal plane position
        eb_pos = pfi.bench.cobras.centers + pfi.bench.cobras.L1 * np.exp(1j * theta)
        fp_pos = eb_pos + pfi.bench.cobras.L2 * np.exp(1j * (theta + phi))

        # Set the broken cobras to home position
        fp_pos[pfi.bench.cobras.hasProblem] = pfi.bench.cobras.home0[pfi.bench.cobras.hasProblem]

        # Pick a few cobras from the good ones and mark them unassigned
        unassigned = np.random.choice(np.arange(pfi.bench.cobras.nCobras)[~pfi.bench.cobras.hasProblem], 10, replace=False)
        fp_pos[unassigned] = TargetGroup.NULL_TARGET_POSITION

        return fp_pos
    
    def test_init(self):
        pfi = self.get_instrument()
        fp_pos = self.generate_targets(pfi)
        targets = TargetGroup(fp_pos, ids=None, priorities=None)
        simulator = CollisionSimulator(pfi, targets)

    def test_calculate_final_fiber_positions(self):
        pfi = self.get_instrument()
        fp_pos = self.generate_targets(pfi)
        targets = TargetGroup(fp_pos, ids=None, priorities=None)
        simulator = CollisionSimulator(pfi, targets)
        simulator._CollisionSimulator__calculate_final_fiber_positions()

    def test_optimize_unassigned_cobra_positions(self):
        pfi = self.get_instrument()
        fp_pos = self.generate_targets(pfi)
        targets = TargetGroup(fp_pos, ids=None, priorities=None)
        simulator = CollisionSimulator(pfi, targets)
        simulator._CollisionSimulator__calculate_final_fiber_positions()
        simulator._CollisionSimulator__optimize_unassigned_cobra_positions()

    def test_calculate_trajectories(self):
        pfi = self.get_instrument()
        fp_pos = self.generate_targets(pfi)
        targets = TargetGroup(fp_pos, ids=None, priorities=None)
        simulator = CollisionSimulator(pfi, targets)
        simulator._CollisionSimulator__calculate_final_fiber_positions()
        simulator._CollisionSimulator__optimize_unassigned_cobra_positions()
        simulator._CollisionSimulator__calculate_trajectories()

    def test_detect_trajectory_collisions(self):
        pfi = self.get_instrument()
        fp_pos = self.generate_targets(pfi)
        targets = TargetGroup(fp_pos, ids=None, priorities=None)
        simulator = CollisionSimulator(pfi, targets)
        simulator._CollisionSimulator__calculate_final_fiber_positions()
        simulator._CollisionSimulator__optimize_unassigned_cobra_positions()
        simulator._CollisionSimulator__calculate_trajectories()
        simulator._CollisionSimulator__detect_trajectory_collisions()

        pass