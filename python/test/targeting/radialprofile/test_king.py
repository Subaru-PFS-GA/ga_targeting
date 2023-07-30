import os
import numpy as np
import numpy.testing as npt

from .test_radialprofile import RadialProfileTest

from pfs.ga.targeting.radialprofile import King
from pfs.ga.targeting.projection import Pointing, WcsProjection, Whitening

class KingTest(RadialProfileTest):
    def get_model(self, obs):
        p = self.get_projection(obs)
        w, mask = self.get_whitening(obs, projection=p)
        xy = w.apply(obs.get_coords(ctype='a')[mask])
        c = xy.mean(axis=0)
        r = King(transformation=w, R_max=6, bins=30, center=c)

        return r

    def test_sample(self):
        r = King(R_max=3)
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
                
    def test_world_to_elliptic(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

        R, phi = r.world_to_elliptic(obs.get_coords())

        self.assertEqual(R.shape, phi.shape)

    def test_elliptic_to_world(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

        ra1, dec1 = obs.get_coords()
        R, phi = r.world_to_elliptic(ra1, dec1)
        ra2, dec2 = r.elliptic_to_world(R, phi, ctype='t')

        npt.assert_almost_equal(ra1, ra2)
        npt.assert_almost_equal(dec1, dec2)

    def test_get_axes_world(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

        a, b, pa = r.get_axes_world()

    def test_get_params_world(self):
        obs = self.load_test_observation()
        r = self.get_model(obs)
        r.histogram(obs)
        r.fit()

        R_c, R_t = r.get_params_world()