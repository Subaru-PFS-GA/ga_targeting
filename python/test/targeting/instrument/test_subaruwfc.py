import os
import numpy as np

from test_base import TestBase

from pfs.ga.common.diagram import FOV
from pfs.ga.common.projection import Pointing
from pfs.ga.targeting.radialprofile import King
from pfs.ga.targeting.instrument import SubaruWFC


class SubaruWfcTest(TestBase):
    def test_world_to_pixel(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        p = Pointing([ ra.mean(), dec.mean() ])
        t = SubaruWFC(p)
        
        (x, y), mask = t.world_to_pixel(ra, dec)
        self.assertIsNotNone(mask)
        self.assertEqual(ra.shape, x.shape)
        self.assertEqual(dec.shape, y.shape)

    def test_world_to_pixel_nans(self):
        pointing = Pointing(0, 0)
        
        radec = np.stack(2 * [np.linspace(-2, 2, 100)], axis=-1)
        
        wfc = SubaruWFC(pointing)
        xy, fov_mask = wfc.world_to_pixel(radec)

    def test_pixel_to_world(self):
        p = Pointing([ 10.0, 10.0 ])
        t = SubaruWFC(p)
        
        [x, y] = np.random.uniform(-200.0, 200.0, size=(2, 1000))

        (ra, dec), fov_mask = t.pixel_to_world(x, y)
        self.assertIsNotNone(fov_mask)
        self.assertEqual(ra.shape, x.shape)
        self.assertEqual(dec.shape, y.shape)

    def test_get_outline_points_pixel(self):
        p = Pointing([ 10.0, 10.0 ])
        t = SubaruWFC(p)

        xy = t._SubaruWFC__get_outline_points_pixel()

    def test_get_outline_points_world(self):
        p = Pointing([ 10.0, 10.0 ])
        t = SubaruWFC(p)

        xy = t._SubaruWFC__get_outline_points_world()

    def test_plot_field_of_view(self):
        obs = self.load_test_observation()
        proj = self.get_projection(obs)
        fov = FOV(projection=proj)
        wfc = SubaruWFC(proj.pointing)

        f, ax = self.get_test_plot(projection=proj.wcs)
        obs.plot(ax, fov)
        wfc.plot_field_of_view(ax, fov)

        self.save_fig(f)

        # TODO: test wrap-around