import os
import numpy as np

from .test_radialprofile import RadialProfileTest

from pfs.ga.targeting.radialprofile import Sersic
from pfs.ga.targeting.projection import Pointing, WcsProjection, Whitening

class SersicTest(RadialProfileTest):
    def get_model(self, obs):
        p = self.get_projection(obs)
        w, mask = self.get_whitening(obs, projection=p)

        r = Sersic(transformation=w, R_max=1)

        return r

    def test_sample(self):
        self.skipTest('Not implemented yet.')

        r = Sersic(R_max=3)
        r.params = (0, 1, 0.5, 1)

        R = r.sample(10)

    def test_histogram(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)

    def test_fit(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

    def test_fit_nomorph(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()
        r.fit_nomorph()

    def test_get_R_halflight(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

        R_c = r.get_R_halflight()
        pass