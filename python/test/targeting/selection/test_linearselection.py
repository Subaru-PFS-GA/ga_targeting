import os
import numpy as np

from test_base import TestBase
from pfs.ga.targeting.photometry import Color
from pfs.ga.targeting.selection import LinearSelection

class LinearSelectionTest(TestBase):
    def test_get_shape(self):
        cmd, photometry = self.get_test_cmd()
        sel = LinearSelection(cmd.axes, np.eye(2), b_min=[17.5, 0.1], b_max=[22.5, 2.0])

        shape = sel.get_shape()

    def test_apply(self):
        cmd, photometry = self.get_test_cmd()
        obs = self.load_test_observation()
        sim = self.load_test_simulation()

        sel = LinearSelection([cmd.axes[1]], 1.0, b_min=17.5, b_max=22.5)
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)

        sel = LinearSelection([cmd.axes[1]], 1.0, b_min=17.5, b_max=22.5)
        mask = sel.apply(sim)
        self.assertEqual(2, mask.ndim)

        sel = LinearSelection([cmd.axes[1]], 1.0, b_min=17.5, b_max=22.5)
        mask = sel.apply(sim, observed=True)
        self.assertEqual(2, mask.ndim)

        sel = LinearSelection([cmd.axes[0]], 1.0, b_min=-1.0, b_max=2.0)
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)

        sel = LinearSelection(cmd.axes, [1, 0.1], -1)
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)