import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing, WcsProjection

class DiagramTest(TestBase):
    def test_plot(self):
        obs = self.load_test_observation()
        f, ax = self.get_test_plot()

        d = self.get_test_diagram()
        ra, dec = obs.get_coords()

        d.plot(ax, ra, dec)

        self.save_fig(f)

    def test_scatter(self):
        obs = self.load_test_observation()
        f, ax = self.get_test_plot()

        d = self.get_test_diagram()
        ra, dec = obs.get_coords()

        d.scatter(ax, ra, dec)

        self.save_fig(f)

    def test_imshow(self):
        self.skipTest()