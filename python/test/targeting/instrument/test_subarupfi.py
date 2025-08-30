import os
import numpy as np
import numpy.testing as npt
import matplotlib.pyplot as plt
from astropy.time import Time

from test_base import TestBase

from ics.cobraCharmer.cobraCoach import engineer
from ics.cobraCharmer.pfiDesign import PFIDesign

from pfs.ga.targeting.config.netflow import NetflowConfig
from pfs.ga.targeting.config.instrument import InstrumentOptionsConfig
from pfs.ga.common.projection import Pointing
from pfs.ga.common.diagram import FOV
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC, CobraAngleFlags

class SubaruPFITest(TestBase):
    def test_init(self):
        pfi = SubaruPFI()

        config = NetflowConfig()
        pfi = SubaruPFI(instrument_options=config.instrument_options)

    def test_load_grand_fiber_map(self):
        pfi = SubaruPFI()
        pfi._SubaruPFI__load_grand_fiber_map()

    def test_get_fiber_map(self):
        pfi = SubaruPFI()
        fiber_map = pfi.fiber_map

    def test_create_configured_bench(self):
        pfi = SubaruPFI()
        pfi._SubaruPFI__create_configured_bench()

    def test_get_bench(self):
        pfi = SubaruPFI()
        bench = pfi.bench

        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)
        bench = pfi.bench

    def test_find_associations(self):
        pfi = SubaruPFI()
        b = pfi.bench
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
        obs.plot(ax, fov)
        wfc.plot_field_of_view(ax, fov)
        pfi.plot_focal_plane(ax, fov, corners=True, c='r', ls='--')
        pfi.plot_cobras(ax, fov, color='b')

        self.save_fig(f)

        # TODO: test wrap-around

    def test_radec_to_altaz(self):
        # TODO: move to telescope class

        pfi = SubaruPFI()
        az, el, inr = pfi.radec_to_altaz(226.3, 67.5, posang=0.0, obs_time=Time("2024-06-10T00:00:00.0Z"))

    def test_get_visibility(self):
        # TODO: move to telescope class

        # Ursa Minor dSph visible
        pfi = SubaruPFI()
        visible, airmass = pfi.get_visibility(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        self.assertTrue(visible)

        # Below horizon
        visible, airmass = pfi.get_visibility(0, 0, posang=0, obs_time=Time("2016-04-03T08:00:00Z"))
        self.assertFalse(visible)

        # Bootes I is visible at a very large airmass
        visible, airmass = pfi.get_visibility(210.025, 14.5, posang=30, obs_time=Time("2025-01-24T11:00:00Z"))
        self.assertTrue(visible)
        self.assertTrue(airmass > 5)

    def test_radec_to_fp_pos(self):
        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)
        pointing = Pointing(226.3, 67.5, posang=0, obs_time=Time("2024-06-10T00:00:00.0Z"))
        fp_pos = pfi.radec_to_fp_pos(np.array(226.3),
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
        pfi = SubaruPFI(instrument_options=instrument_options)

        s = np.s_[:]
        cidx = np.arange(pfi.bench.cobras.nCobras)[s]
        centers = pfi.bench.cobras.centers[s][:, None]

        # Generate random focal plane positions around each cobra
        # Radius is very slightly beyond the maximum reach of the cobra
        fp_pos = self.generate_random_fp_pos(centers)

        theta, phi, d, eb_pos, flags = pfi.fp_pos_to_cobra_angles(fp_pos, cidx)

        self.assertEqual(theta.shape, (2394, 100, 2))
        self.assertEqual(phi.shape, (2394, 100, 2))

        # Make sure there's only secondary solution of there is a first one
        self.assertTrue(((flags[..., 1] & CobraAngleFlags.SOLUTION_OK) & ~(flags[..., 0] & CobraAngleFlags.SOLUTION_OK)).sum() == 0)

        # Make sure the two solutions evaluate to the same position
        fp_pos2 = pfi.cobra_angles_to_fp_pos(theta, phi, cidx)

        mask = (flags & CobraAngleFlags.SOLUTION_OK != 0)
        npt.assert_almost_equal(fp_pos[mask[..., 0]], fp_pos2[..., 0][mask[..., 0]])
        npt.assert_almost_equal(fp_pos[mask[..., 1]], fp_pos2[..., 1][mask[..., 1]])

        # Compare to CobraOps implementation

        # First solutions are the same for non-broken cobras
        # Many second solutions are NaN, so not testing them for now
        for j in range(fp_pos.shape[1]):
            theta2, phi2, flags2 = pfi.cobra_coach.pfi.positionsToAngles(pfi.cobra_coach.allCobras[s], fp_pos[:, j])
            for i in range(2):
                npt.assert_almost_equal(theta[:, j, i][pfi.bench.cobras.isGood[s]], theta2[:, i][pfi.bench.cobras.isGood[s]])
                npt.assert_almost_equal(phi[:, j, i][pfi.bench.cobras.isGood[s]], phi2[:, i][pfi.bench.cobras.isGood[s]])
                
                # Flags don't match by definition
                # self.assertTrue((flags[:, j, i] == flags2[:, i]).all())

        # Verify that the elbow positions are correct
        L2 = pfi.bench.cobras.L2[cidx][..., None]
        theta0 = pfi.bench.cobras.tht0[cidx][..., None]

        solution_ok = (flags[..., 0] & CobraAngleFlags.SOLUTION_OK) != 0
        npt.assert_almost_equal((np.abs(eb_pos[..., 0] - fp_pos) - L2)[solution_ok], 0)

        solution_ok = (flags[..., 1] & CobraAngleFlags.SOLUTION_OK) != 0
        npt.assert_almost_equal((np.abs(eb_pos[..., 1] - fp_pos) - L2)[solution_ok], 0)

    def test_cobra_angles_to_fp_pos(self):
        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)

        s = np.s_[:120]
        cidx = np.arange(pfi.bench.cobras.nCobras)[s]
        centers = pfi.bench.cobras.centers[s][:, None]

        fp_pos = self.generate_random_fp_pos(centers)

        theta, phi, d, eb_pos, flags = pfi.fp_pos_to_cobra_angles(fp_pos, cidx)
        fp_pos2 = pfi.cobra_angles_to_fp_pos(theta[..., 0], phi[..., 0], cidx)

        pass

    def test_simulate_trajectories(self):
        # This is my version

        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)

        cidx = np.where(pfi.calib_model.status == PFIDesign.COBRA_OK_MASK)[0][:120]
        centers = pfi.bench.cobras.centers[cidx][:, None]

        fp_pos = self.generate_random_fp_pos(centers, n=10)
        cobraidx = np.broadcast_to(cidx[..., None], fp_pos.shape)
        theta, phi, d, eb_pos, flags = pfi.fp_pos_to_cobra_angles(fp_pos, cidx)

        # Send all cobras to home position
        cobra_state = pfi.create_cobra_state(cobraidx)
        pfi.set_cobra_state_to_home(cobra_state)
        
        # Simulate the trajectories from the current state to the target state
        # Use the first solution for the angles
        pfi.simulate_trajectories(cobra_state, theta[..., 0], phi[..., 0],
                                  use_scaling=pfi.cobra_coach.useScaling, 
                                  max_segments=pfi.cobra_coach.maxSegments,
                                  max_total_steps=pfi.cobra_coach.maxTotalSteps)

    def test_simulate_trajectories_cobraCharmer(self):
        # This is how it's done in the original code


        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)
        cc = pfi.cobra_coach

        s = np.s_[:120]
        cidx = np.arange(pfi.bench.cobras.nCobras)[s]
        centers = pfi.bench.cobras.centers[s][:, None]

        fp_pos = self.generate_random_fp_pos(centers)

        # Calculate the final theta and phi angles for the good cobras
        timeStep=20
        maxSteps = 2000
        nCobras = cc.nCobras
        goodCobras = np.full(nCobras, False)
        goodCobras[cc.goodIdx] = True

        thetaAngles, phiAngles, _ = cc.pfi.positionsToAngles(cc.allCobras[s][goodCobras[s]], fp_pos[goodCobras[s], 0])

        # Select the first angles solution
        thetaAngles = thetaAngles[:, 0]
        phiAngles = phiAngles[:, 0]

        # Initialize the engineer module
        engineer.setCobraCoach(cc)
        engineer.setConstantOntimeMode(maxSteps=maxSteps)

        # TODO: coordinates must be local!
        #       compare with CollisionSimulator2
        # raise NotImplementedError()

        # Calculate the cobra trajectories
        trajectories, _ = engineer.createTrajectory(
            np.where(goodCobras[s])[0], thetaAngles, phiAngles,
            tries=8, twoSteps=True, threshold=20.0, timeStep=timeStep)

        # Calculate the fiber and elbow positions along the cobra trajectories
        fiberPositions = trajectories.calculateFiberPositions(cc)
        elbowPositions = trajectories.calculateElbowPositions(cc)
        
        # self.nSteps = self.fiberPositions.shape[1]

        pass

    def test_batch_interpolate(self):
        x = np.random.uniform(0, 113, (120, 80, 100))
        xp = np.broadcast_to(np.arange(113), (120, 80, 113))
        fp = np.zeros((120, 80, 113))

        yp = SubaruPFI._SubaruPFI__batch_interp(None, x, xp, fp)
        self.assertEqual(yp.shape, (120, 80, 100))

    def test_home_position(self):
        # Verify how the home position is calculated from the angles

        instrument_options = InstrumentOptionsConfig.from_dict({'layout': 'calibration'})
        pfi = SubaruPFI(instrument_options=instrument_options)

        # pfi.bench.cobras.home0
        # pfi.bench.cobras.home1
        # pfi.bench.cobras.L1
        # pfi.bench.cobras.L2
        # pfi.bench.cobras.tht0
        # pfi.bench.cobras.tht1
        # pfi.bench.cobras.phiIn
        # pfi.bench.cobras.phiOut
        # pfi.bench.cobras.phiHome

        # phi.calib_model.centers
        # phi.calib_model.L1
        # phi.calib_model.L2
        # phi.calib_model.phiIn
        # phi.calib_model.phiOut
        # phi.calib_model.tht0
        # phi.calib_model.tht1

        npt.assert_equal(pfi.bench.cobras.centers, pfi.calib_model.centers)
        npt.assert_equal(pfi.bench.cobras.L1, pfi.calib_model.L1)
        npt.assert_equal(pfi.bench.cobras.L2, pfi.calib_model.L2)
        npt.assert_equal(pfi.bench.cobras.tht0, pfi.calib_model.tht0)
        npt.assert_equal(pfi.bench.cobras.tht1, pfi.calib_model.tht1)
        npt.assert_equal(pfi.bench.cobras.phiIn, pfi.calib_model.phiIn)
        npt.assert_equal(pfi.bench.cobras.phiOut, pfi.calib_model.phiOut)
        
        # CobraOps version
        # The home position angles are given in global angles in cobraOps
        phi_home = np.maximum(-np.pi, pfi.bench.cobras.phiIn) + 0.00001
        home0_cobraops = pfi.bench.cobras.centers + \
                         pfi.bench.cobras.L1 * np.exp(1j * pfi.bench.cobras.tht0) + \
                         pfi.bench.cobras.L2 * np.exp(1j * (pfi.bench.cobras.tht0 + phi_home))
        home1_cobraops = pfi.bench.cobras.centers + \
                         pfi.bench.cobras.L1 * np.exp(1j * pfi.bench.cobras.tht1) + \
                         pfi.bench.cobras.L2 * np.exp(1j * (pfi.bench.cobras.tht1 + phi_home))

        # home0 and home1 are removed from cobraOps
        # npt.assert_almost_equal(home0_cobraops, pfi.bench.cobras.home0)
        # npt.assert_almost_equal(home1_cobraops[~pfi.bench.cobras.hasProblem], pfi.bench.cobras.home1[~pfi.bench.cobras.hasProblem])
        # npt.assert_almost_equal(home0_cobraops[pfi.bench.cobras.hasProblem], pfi.bench.cobras.home1[pfi.bench.cobras.hasProblem])

        # CobraCoach version (trajectory mode, this is what the simulator uses)
        # There are two modes for theta, CW and CWW
        # These are implemented in CobraCoach.moveToHome
        theta_home_CWW = np.zeros_like(pfi.calib_model.centers)
        theta_home_CW = (pfi.calib_model.tht1 - pfi.calib_model.tht0 + np.pi) % (2 * np.pi) + np.pi
        phi_home = np.zeros_like(pfi.calib_model.centers)
        home_CWW = pfi.calib_model.centers + \
                   pfi.calib_model.L1 * np.exp(1j * theta_home_CWW) + \
                   pfi.calib_model.L2 * np.exp(1j * (phi_home + pfi.calib_model.tht0))
        home_CW = pfi.calib_model.centers + \
                  pfi.calib_model.L1 * np.exp(1j * theta_home_CW) + \
                  pfi.calib_model.L2 * np.exp(1j * (phi_home + pfi.calib_model.tht0)) 

        #

        pass