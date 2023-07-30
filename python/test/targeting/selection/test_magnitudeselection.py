import os
import numpy as np

from test_base import TestBase
from pfs.ga.targeting.photometry import Color
from pfs.ga.targeting.selection import MagnitudeSelection

class MagnitudeSelectionTest(TestBase):
    def test_apply(self):
        cmd, photometry = self.get_test_cmd()
        obs = self.load_test_observation()
        sim = self.load_test_simulation()

        sel = MagnitudeSelection(cmd.axes[1], 18.5, 22.5)
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)

        sel = MagnitudeSelection(cmd.axes[1], 17, 22.5)
        mask = sel.apply(sim)
        self.assertEqual(2, mask.ndim)