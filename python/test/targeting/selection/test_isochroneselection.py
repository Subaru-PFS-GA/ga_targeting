import os
import numpy as np


from pfs.ga.common.diagram import CMD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Color, Magnitude, Photometry, photometry
from pfs.ga.targeting import Isochrone
from pfs.ga.targeting.selection import IsochroneSelection

from test_base import TestBase

class IsochroneSelectionTest(TestBase):
    def get_selection(self, cmd, iso):
        sel = IsochroneSelection(iso, 
            cmd.axes,
            selection_axis=0, selection_direction='+',
            DM=19.2, error_sigma=[-2, 0])

        return sel

    def test_apply(self):
        obs = self.load_test_observation()
        cmd, photometry = self.get_test_cmd()
        iso = self.get_test_isochrone()
        sel = self.get_selection(cmd, iso)
        
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)
