import os
import numpy as np
import matplotlib.pyplot as plt

from ics.cobraOps.Bench import Bench

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC
from pfs.ga.targeting.diagram import FOV

class SubaruPFITest(TestBase):
    def test_find_associations(self):
        b = Bench()
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        mask = np.full_like(ra, True, dtype=bool)

        wfc = SubaruWFC(pointing=Pointing(227.2420, 67.2221))
        pfi = SubaruPFI(projection=wfc)
        assoc = pfi.find_associations(ra, dec, mask=mask)

        self.assertEqual(len(b.cobras.centers), len(assoc.assoc))

    def test_plot_field_of_view(self):
        obs = self.load_test_observation()
        proj = self.get_projection(obs)
        fov = FOV(projection=proj)
        wfc = SubaruWFC(proj.pointing)
        pfi = SubaruPFI(wfc)

        f, ax = self.get_test_plot(projection=proj.wcs)
        fov.plot_observation(ax, obs)
        wfc.plot_field_of_view(ax, fov)
        pfi.plot_focal_plane(ax, fov, corners=True, c='r', ls='--')
        pfi.plot_cobras(ax, fov, color='b')

        self.save_fig(f)

        # TODO: test wrap-around