import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.common.projection import Pointing
from pfs.ga.common.diagram import FP
from pfs.ga.targeting.instrument import SubaruWFC, SubaruPFI

class FPTest(TestBase):
    def get_projection(self, obs):
        ra, dec = obs.get_coords()
        p = SubaruWFC(Pointing(ra.mean(), dec.mean()))
        return p

    def test_plot_world(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        wfc = self.get_projection(obs)
        fp = FP(projection=wfc)
        
        f, ax = self.get_test_plot()
        fp.plot(ax, ra, dec, native_frame='world', fmt='.')
        wfc.plot_focal_plane(ax, fp)
        ax.grid()
        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        
        self.save_fig(f)

    def test_plot_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fp = FP(projection=p)
        
        f, ax = self.get_test_plot()
        fp.plot(ax, ra, dec, native_frame='pixel')
        ax.grid()
        
        self.save_fig(f)

    def test_scatter_world(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fp = FP(projection=p)
        
        f, ax = self.get_test_plot()
        fp.scatter(ax, ra, dec, native_frame='world', s=1)
        ax.set_xlim(-250, 250)
        ax.set_ylim(-250, 250)
        ax.grid()
        
        self.save_fig(f)

    def test_scatter_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fp = FP(projection=p)
        
        f, ax = self.get_test_plot()
        fp.scatter(ax, ra, dec, native_frame='pixel')
        ax.grid()
        
        self.save_fig(f)

    def test_plot_observation(self):
        obs = self.load_test_observation()
        p = self.get_projection(obs)
        fp = FP(projection=p)

        f, ax = self.get_test_plot()
        obs.plot(ax, fp)

        self.save_fig(f)

    def test_plot_instrument(self):
        obs = self.load_test_observation()
        wfc = self.get_projection(obs)
        pfi = SubaruPFI(projection=wfc)
        fp = FP(projection=wfc)
        
        f, ax = self.get_test_plot()
        wfc.plot_focal_plane(ax, fp)
        pfi.plot_focal_plane(ax, fp, corners=True)
        pfi.plot_focal_plane(ax, fp, blocks=True, c='blue')
        pfi.plot_cobras(ax, fp)
        
        ax.grid()
        ax.set_xlim(-260, 260)
        ax.set_ylim(-260, 260)
        
        self.save_fig(f)