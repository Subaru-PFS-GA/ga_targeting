import os
import numpy as np
import numpy.testing as npt
import matplotlib.pyplot as plt
from astropy.time import Time

from test_base import TestBase

from pfs.ga.targeting.config import NetflowConfig
from pfs.ga.targeting.config.instrumentoptionsconfig import InstrumentOptionsConfig
from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC, CobraAngleFlags
from pfs.ga.targeting.diagram import FOV

class SubaruPFITest(TestBase):
    def test_init(self):
        inst = SubaruPFI()

        config = NetflowConfig()
        inst = SubaruPFI(instrument_options=config.instrument_options)

    def test_load_grand_fiber_map(self):
        inst = SubaruPFI()
        inst._SubaruPFI__load_grand_fiber_map()

    def test_get_fiber_map(self):
        inst = SubaruPFI()
        fiber_map = inst.fiber_map

    def test_create_default_bench(self):
        inst = SubaruPFI()
        inst._SubaruPFI__create_default_bench()

    def test_create_configured_bench(self):
        inst = SubaruPFI()
        inst._SubaruPFI__create_configured_bench()

    def test_get_bench(self):
        inst = SubaruPFI()
        bench = inst.bench

        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'full'})
        inst = SubaruPFI(instrument_options=instrument_options)
        bench = inst.bench

        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        inst = SubaruPFI(instrument_options=instrument_options)
        bench = inst.bench

    def test_find_associations(self):
        inst = SubaruPFI()
        b = inst.bench
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

    def test_radec_to_altaz(self):
        # TODO: move to telescope class

        inst = SubaruPFI()
        az, el, inr = inst.radec_to_altaz(226.3, 67.5, posang=0.0, obs_time=Time("2024-06-10T00:00:00.0Z"))

    def test_get_visibility(self):
        # TODO: move to telescope class

        # Ursa Minor dSph visible
        inst = SubaruPFI()
        visible, airmass = inst.get_visibility(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        self.assertTrue(visible)

        # Below horizon
        visible, airmass = inst.get_visibility(0, 0, posang=0, obs_time=Time("2016-04-03T08:00:00Z"))
        self.assertFalse(visible)

        # Bootes I is visible at a very large airmass
        visible, airmass = inst.get_visibility(210.025, 14.5, posang=30, obs_time=Time("2025-01-24T11:00:00Z"))
        self.assertTrue(visible)
        self.assertTrue(airmass > 5)

    def test_radec_to_fp_pos(self):
        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        inst = SubaruPFI(instrument_options=instrument_options)
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        fp_pos = inst.radec_to_fp_pos(np.array(226.3),
                                      np.array(67.5),
                                      pointing=pointing)
        
    def generate_random_fp_pos(self, centers, n=100):
        np.random.seed(0)
        phi = np.random.uniform(0, 2*np.pi, size=n)
        r = np.sqrt(np.random.uniform(0, 1.0, size=n)) * 5.6
        fp_pos = centers + r * np.cos(phi) + r * np.sin(phi) * 1j
        return fp_pos
        
    def test_fp_pos_to_cobra_angles(self):
        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        inst = SubaruPFI(instrument_options=instrument_options)

        s = np.s_[:]
        cidx = np.arange(inst.bench.cobras.nCobras)[s]
        centers = inst.bench.cobras.centers[s][:, None]

        # Generate random focal plane positions around each cobra
        # Radius is very slightly beyond the maximum reach of the cobra
        fp_pos = self.generate_random_fp_pos(centers)

        theta, phi, d, eb_pos, flags = inst.fp_pos_to_cobra_angles(fp_pos, cidx)

        self.assertEqual(theta.shape, (2394, 100, 2))
        self.assertEqual(phi.shape, (2394, 100, 2))

        # Make sure there's only secondary solution of there is a first one
        self.assertTrue(((flags[..., 1] & CobraAngleFlags.SOLUTION_OK) & ~(flags[..., 0] & CobraAngleFlags.SOLUTION_OK)).sum() == 0)

        # Make sure the two solutions evaluate to the same position
        fp_pos2 = inst.cobra_angles_to_fp_pos(theta, phi, cidx)

        mask = (flags & CobraAngleFlags.SOLUTION_OK != 0)
        npt.assert_almost_equal(fp_pos[mask[..., 0]], fp_pos2[..., 0][mask[..., 0]])
        npt.assert_almost_equal(fp_pos[mask[..., 1]], fp_pos2[..., 1][mask[..., 1]])

        # Compare to CobraOps implementation

        # First solutions are the same for non-broken cobras
        # Many second solutions are NaN, so not testing them for now
        for j in range(fp_pos.shape[1]):
            theta2, phi2, flags2 = inst.cobra_coach.pfi.positionsToAngles(inst.cobra_coach.allCobras[s], fp_pos[:, j])
            for i in range(2):
                npt.assert_almost_equal(theta[:, j, i][~inst.bench.cobras.hasProblem[s]], theta2[:, i][~inst.bench.cobras.hasProblem[s]])
                npt.assert_almost_equal(phi[:, j, i][~inst.bench.cobras.hasProblem[s]], phi2[:, i][~inst.bench.cobras.hasProblem[s]])
                
                # Flags don't match by definition
                # self.assertTrue((flags[:, j, i] == flags2[:, i]).all())

        # Verify that the elbow positions are correct
        L2 = inst.bench.cobras.L2[cidx][..., None]
        theta0 = inst.bench.cobras.tht0[cidx][..., None]

        solution_ok = (flags[..., 0] & CobraAngleFlags.SOLUTION_OK) != 0
        npt.assert_almost_equal((np.abs(eb_pos[..., 0] - fp_pos) - L2)[solution_ok], 0)

        solution_ok = (flags[..., 1] & CobraAngleFlags.SOLUTION_OK) != 0
        npt.assert_almost_equal((np.abs(eb_pos[..., 1] - fp_pos) - L2)[solution_ok], 0)

    def test_cobra_angles_to_fp_pos(self):
        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        inst = SubaruPFI(instrument_options=instrument_options)

        s = np.s_[:120]
        cidx = np.arange(inst.bench.cobras.nCobras)[s]
        centers = inst.bench.cobras.centers[s][:, None]

        fp_pos = self.generate_random_fp_pos(centers)

        theta, phi, d, eb_pos, flags = inst.fp_pos_to_cobra_angles(fp_pos, cidx)
        fp_pos2 = inst.cobra_angles_to_fp_pos(theta[..., 0], phi[..., 0], cidx)

        pass

    # def test_verify_cobra_angles(self):
    #     instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
    #     inst = SubaruPFI(instrument_options=instrument_options)

    #     s = np.s_[:120]
    #     cidx = np.arange(inst.bench.cobras.nCobras)[s]
    #     centers = inst.bench.cobras.centers[s][:, None]

    #     # Generate random focal plane positions around each cobra
    #     # Radius is very slightly beyond the maximum reach of the cobra

    #     phi = np.random.uniform(0, 2*np.pi, size=(1, 100))
    #     r = np.sqrt(np.random.uniform(0, 1.0, size=(1, 100))) * 5.6
    #     fp_pos = centers + r * np.cos(phi) + r * np.sin(phi) * 1j

    #     theta, phi = inst.fp_pos_to_cobra_angles(fp_pos, cidx)
    #     mask = inst.verify_cobra_angles(theta, phi, cidx)

    #     pass

