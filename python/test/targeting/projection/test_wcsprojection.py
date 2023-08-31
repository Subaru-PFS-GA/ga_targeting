import os
import numpy as np

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing, WcsProjection

class WcsProjectionTest(TestBase):
    def test_world_to_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        
        p = WcsProjection(Pointing(ra.mean(), dec.mean()), proj='TAN')
        
        xy, mask = p.world_to_pixel(obs)
        self.assertEqual((ra.shape[0], 2), xy.shape)

        (x, y), mask = p.world_to_pixel(obs.get_coords())
        self.assertEqual(ra.shape, x.shape)
        self.assertEqual(dec.shape, y.shape)

    def test_pixel_to_world(self):
        xy = np.random.uniform(-1, 1, size=(1000, 2))
        x = np.random.uniform(-1, 1, size=(1000,))
        y = np.random.uniform(-1, 1, size=(1000,))
        
        p = WcsProjection(Pointing(10, 20), proj='TAN')
        radec, mask = p.pixel_to_world(xy)
        self.assertEqual(xy.shape, radec.shape)

        p = WcsProjection(Pointing(10, 20), proj='TAN')
        (ra, dec), mask = p.pixel_to_world(x, y)
        self.assertEqual(x.shape, ra.shape)
        self.assertEqual(y.shape, dec.shape)
