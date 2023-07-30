import os
import numpy as np

from test_base import TestBase

from pfs.ga.targeting.radialprofile.radialprofile import RadialProfile
from pfs.ga.targeting.projection import Pointing, WcsProjection, Whitening

class RadialProfileTest(TestBase):
    def get_projection(self, obs):
        ra, dec = obs.get_coords()
        p = WcsProjection(Pointing(ra.mean(), dec.mean()), proj='TAN')
        return p

    def get_whitening(self, obs, projection=None):
        w = Whitening(projection=projection)
        mask = w.create(obs)
        return w, mask