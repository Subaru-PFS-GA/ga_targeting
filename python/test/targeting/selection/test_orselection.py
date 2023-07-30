import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.targeting.selection import OrSelection

class OrSelectionTest(TestBase):
    def test_apply(self):
        cmd, photometry = self.get_test_cmd()
        obs = self.load_test_observation()

        sel = OrSelection([
            self.get_test_magnitude_selection(cmd),
            self.get_test_color_selection(cmd),
        ])

        mask = sel.apply(obs)
        
