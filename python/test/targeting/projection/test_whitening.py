import os
import numpy as np

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing, WcsProjection, Whitening

class WcsProjectionTest(TestBase):
    def get_projection(self, obs):
        ra, dec = obs.get_coords()
        p = WcsProjection(Pointing(ra.mean(), dec.mean()), proj='TAN')
        return p

    def get_whitening(self, obs, projection=None):
        w = Whitening(projection=projection)
        w.create(obs)
        return w

    def test_create_noproj(self):
        obs = self.load_test_observation()
        w = self.get_whitening(obs)

    def test_create_proj(self):
        obs = self.load_test_observation()
        p = self.get_projection(obs)
        w = self.get_whitening(obs, projection=p)

    def test_apply_noproj(self):
        obs = self.load_test_observation()
        w = self.get_whitening(obs)
        
        wxy = w.apply(obs)

    def test_apply_proj(self):
        obs = self.load_test_observation()
        p = self.get_projection(obs)
        w = self.get_whitening(obs, projection=p)
        
        wxy = w.apply(obs)

    def test_reverse_noproj(self):
        obs = self.load_test_observation()
        w = self.get_whitening(obs)
        wxy = w.apply(obs)
        
        xy = w.reverse(wxy)

    def test_reverse_proj(self):
        obs = self.load_test_observation()
        p = self.get_projection(obs)
        w = self.get_whitening(obs, projection=p)
        wxy = w.apply(obs)
        
        xy = w.reverse(wxy)