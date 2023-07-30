import os
import numpy as np
from shapely.geometry import Polygon

from test_base import TestBase
from pfs.ga.targeting.photometry import Color
from pfs.ga.targeting.selection import PolygonSelection

class PolygonSelectionTest(TestBase):
    def get_test_points(self):
        return [[0.1, 16.1],
                [2.0, 16.1],
                [2.0, 23.5],
                [0.1, 23.5],
                [0.1, 16.1]]

    def test_apply_observation(self):
        cmd, photometry = self.get_test_cmd()
        obs = self.load_test_observation()
        points = self.get_test_points()

        sel = PolygonSelection(cmd.axes, points)
        mask = sel.apply(obs)
        self.assertEqual(1, mask.ndim)

    def test_apply_simulation(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()
        points = self.get_test_points()
    
        sel = PolygonSelection(cmd.axes, points)
        mask = sel.apply(sim, mask=np.s_[:100, :])
        self.assertEqual(2, mask.ndim)

        sel = PolygonSelection(cmd.axes, points)
        mask = sel.apply(sim, observed=True, mask=np.s_[:100, :])
        self.assertEqual(2, mask.ndim)