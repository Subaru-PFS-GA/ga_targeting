import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC
from pfs.ga.targeting.diagram import FOV

class SubaruPFITest(TestBase):
    def test_init(self):
        inst = SubaruPFI()
        inst = SubaruPFI(instrument_options={})

    def test_get_grand_fiber_map(self):
        inst = SubaruPFI()
        inst._SubaruPFI__get_grand_fiber_map()

    def test_get_fiber_map(self):
        inst = SubaruPFI()
        inst.get_fiber_map()

    def test_get_configured_bench(self):
        inst = SubaruPFI()
        inst._Instrument__get_configured_bench()

    def test_get_bench(self):
        inst = SubaruPFI()
        inst.get_bench()

        inst = SubaruPFI(instrument_options={'layout': 'full'})
        inst.get_bench()

        inst = SubaruPFI(instrument_options={'layout': 'calibration'})
        inst.get_bench()

    def test_radec_to_altaz(self):
        inst = SubaruPFI()
        inst.radec_to_altaz(226.3, 67.5, posang=0.0, obs_time=Time("2016-04-03T08:00:00Z"))

    def test_radec_to_fp_pos(self):
        inst = SubaruPFI()
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        inst.radec_to_fp_pos(pointing,
                             np.array(226.3),
                             np.array(67.5))

    def test_find_associations(self):
        inst = SubaruPFI()
        b = int.bench
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

    def test_nf_get_closest_dots(self):
        pfi = SubaruPFI()
        dots = pfi.nf_get_closest_dots()
        self.assertEqual(2394, len(dots))