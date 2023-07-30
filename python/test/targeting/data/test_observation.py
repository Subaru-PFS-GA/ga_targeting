import os
import numpy as np

from test_base import TestBase

class ObservationTest(TestBase):
    def test_get_magnitude(self):
        cmd, _ = self.get_test_cmd()
        obs = self.load_test_observation()

        m, s = obs.get_magnitude(cmd.axes[1].magnitude)
        self.assertIsNotNone(m)
        self.assertIsNotNone(s)

        m, s = obs.get_magnitude(cmd.axes[1].magnitude, dered=False)
        self.assertIsNotNone(m)
        self.assertIsNotNone(s)