import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.targeting.diagram import FOV

class FOVTest(TestBase):
    def test_plot(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.plot(ax, ra, dec)
        ax.grid()

        self.save_fig(f)

    def test_plot_world(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.plot(ax, ra, dec, native_frame='world')
        ax.grid()
        
        self.save_fig(f)

    def test_plot_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.plot(ax, ra, dec, native_frame='pixel')
        ax.grid()
        
        self.save_fig(f)

    def test_scatter(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.scatter(ax, ra, dec)
        ax.grid()
        
        self.save_fig(f)

    def test_scatter_world(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.scatter(ax, ra, dec, native_frame='world')
        ax.grid()
        
        self.save_fig(f)

    def test_scatter_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = self.get_projection(obs)
        fov = FOV(projection=p)
        
        f, ax = self.get_test_plot(projection=p.wcs)
        fov.scatter(ax, ra, dec, native_frame='pixel')
        ax.grid()
        
        self.save_fig(f)

    def test_plot_observation(self):
        obs = self.load_test_observation()
        p = self.get_projection(obs)
        fov = FOV(projection=p)

        f, ax = self.get_test_plot(projection=p.wcs)
        fov.plot_observation(ax, obs)

        self.save_fig(f)