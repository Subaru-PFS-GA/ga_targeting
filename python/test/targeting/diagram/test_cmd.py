import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing, WcsProjection
from pfs.ga.targeting.diagram import CMD

class CMDTest(TestBase):
    def test_plot_selection(self):
        raise NotImplementedError()

    def test_plot_isochrone(self):
        iso = self.get_test_isochrone()
        f, ax = self.get_test_plot()

        cmd, _ = self.get_test_cmd()
        cmd.plot_isochrone(ax, iso)

        self.save_fig(f)

    def test_plot_catalog(self):
        obs = self.load_test_observation()
        f, ax = self.get_test_plot()

        cmd, _ = self.get_test_cmd()
        cmd.plot_catalog(ax, obs)

        self.save_fig(f)

    def test_plot_observation(self):
        obs = self.load_test_observation()
        f, ax = self.get_test_plot()

        cmd, _ = self.get_test_cmd()
        cmd.plot_observation(ax, obs)

        self.save_fig(f)

    def test_plot_simulation(self):
        sim = self.load_test_simulation()
        f, ax = self.get_test_plot()

        cmd, _ = self.get_test_cmd()
        cmd.plot_simulation(ax, sim, apply_categories=True, s=np.s_[::30])

        self.save_fig(f)

    def test_plot_probability_map(self):
        raise NotImplementedError()