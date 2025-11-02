import os
from collections import defaultdict
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.cm import get_cmap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from sklearn.neighbors import KDTree
from scipy.spatial import distance_matrix

import pfs.utils
from pfs.utils.coordinates import DistortionCoefficients as DCoeff
from pfs.utils.coordinates.CoordTransp import CoordinateTransform
import pfs.instdata
from ics.cobraOps.Bench import Bench
from pfs.utils.fiberids import FiberIds
from pfs.utils.butler import Butler

from pfs.ga.common.util import *
from pfs.ga.common.util.config import *
from pfs.ga.common.diagram import styles
from pfs.ga.common.projection import Pointing

from ..external import BlackDotsCalibrationProduct
from ..external import PFIDesign
from ..external import CobraCoach
from ..allocation.associations import Associations
from ..allocation.fiberallocator import FiberAllocator
from .instrument import Instrument
from .subaruwfc import SubaruWFC
from .cobraangleflags import CobraAngleFlags

from ..setup_logger import logger

class SubaruPFI(Instrument, FiberAllocator):
    """
    Implements a wrapper around the Subaru PFI instrument configuration.

    Configuration steps are adopted from the PFS design tool:
    ets_pointing/src/pfs_design_tool/pointing_utils/nfutils.py

    Cobra positions video: https://www.newscaletech.com/video-cobra-fiber-positioner-for-nasa-jpl/
    """

    # Path to grand fiber map file
    DEFAULT_FIBERIDS_PATH = os.path.join(os.path.dirname(pfs.utils.__file__), '../../../data/fiberids')
    
    # Path to instrument configuration
    # PFS_INSTDATA_DIR is set to this path
    DEFAULT_INSTDATA_PATH = os.path.join(os.path.dirname(pfs.instdata.__file__), '../../../')

    # Path to black dots calibration file
    DEFAULT_BLACK_DOTS_PATH = os.path.join(DEFAULT_INSTDATA_PATH, "data/pfi/dot", "black_dots_mm.csv")

    # Path to the cobra coach directory
    DEFAULT_COBRA_COACH_DIR = os.path.expandvars('/tmp/$USER/cobraCoach')

    # IDs of cobras in corners
    CORNERS = [[742, 796, 56, 2338, 2392, 1652, 1540, 1594, 854]]

    # IDs of cobras in cobra block corners
    BLOCKS = [[742, 796, 56, 2],
                [1598, 2338, 2392, 1652],
                [854, 800, 1540, 1594]]

    fiber_map_cache = {}
    bench_cache = {}

    def __init__(self, projection=None, instrument_options=None, use_cached_bench=True, orig=None):
        from ..config.instrument.instrumentoptionsconfig import InstrumentOptionsConfig
        
        Instrument.__init__(self, orig=orig)
        FiberAllocator.__init__(self, orig=orig)

        if isinstance(instrument_options, dict):
            instrument_options = InstrumentOptionsConfig.from_dict(instrument_options)

        if not isinstance(orig, SubaruPFI):
            self.__projection = projection or SubaruWFC(Pointing(0, 0))
            self.__instrument_options = instrument_options if instrument_options is not None else InstrumentOptionsConfig.default()
        else:
            self.__projection = projection or orig.__projection
            self.__instrument_options = orig.__instrument_options

        self.__fiber_map, self.__blocked_fibers = self.__load_grand_fiber_map()
        self.__bench, self.__cobra_coach, self.__calib_model = self.__create_bench(use_cached_bench=use_cached_bench)

    #region Properties

    def __get_fiber_map(self):
        return self.__fiber_map
    
    fiber_map = property(__get_fiber_map)

    def __get_blocked_fibers(self):
        return self.__blocked_fibers
    
    blocked_fibers = property(__get_blocked_fibers)

    def __get_calib_model(self):
        return self.__calib_model
    
    calib_model = property(__get_calib_model)

    def __get_bench(self):
        return self.__bench

    bench = property(__get_bench)

    def __get_cobra_coach(self):
        return self.__cobra_coach
    
    cobra_coach = property(__get_cobra_coach)

    def __get_cobra_count(self):
        return self.__bench.cobras.nCobras

    cobra_count = property(__get_cobra_count)    

    #endregion
        
    def __get_instrument_option(self, option, default=None):
        """
        Return an option from the instrument_options dictionary, if it exists,
        otherwise return the default value.
        """

        if option is not None:
            return option
        else:
            return default
        
    def __load_grand_fiber_map(self):
        # Load the grand fiber map
        fiberids_path = self.__get_instrument_option(self.__instrument_options.fiberids_path, SubaruPFI.DEFAULT_FIBERIDS_PATH)

        if fiberids_path in SubaruPFI.fiber_map_cache:
            logger.info(f"Getting the grand fiber map from cache.")

            (fiber_map, blocked_fibers) = SubaruPFI.fiber_map_cache[fiberids_path]
        else:
            fiber_map = FiberIds(path=fiberids_path)
            
            logger.info(f"Loading the grand fiber map from `{os.path.abspath(fiberids_path)}`")

            # Load the list of blocked fibers
            config_root = os.path.join(self.__get_instrument_option(self.__instrument_options.instdata_path, SubaruPFI.DEFAULT_INSTDATA_PATH), 'data')
            butler = Butler(configRoot=config_root)
            blocked_fibers = butler.get('fiberBlocked').set_index('fiberId')

            logger.info(f"Loaded configuration for {len(fiber_map.fiberId)} fibers.")
            logger.info(f"Number of fibers connected to cobras: {(fiber_map.cobraId != fiber_map.MISSING_VALUE).sum()}")
            logger.info(f"Number of science fibers: {(fiber_map.scienceFiberId < fiber_map.ENGINEERING).sum()}")

            SubaruPFI.fiber_map_cache[fiberids_path] = (fiber_map, blocked_fibers)

        return fiber_map, blocked_fibers

    def __load_calib_model(self):
        """
        Load the instrument calibration model directly.

        We don't use this but load the calibration model through cobra coach instead.
        """
        config_root = os.path.join(self.__get_instrument_option(self.__instrument_options.instdata_path, SubaruPFI.DEFAULT_INSTDATA_PATH), 'data')
        butler = Butler(configRoot=config_root)

        # This call might trigger a few division by zero warnings but they are supposed to be OK
        calib_model = butler.get('moduleXml', moduleName='ALL', version='')
        return calib_model
    
    def __create_bench(self, use_cached_bench=True):
        layout = self.__get_instrument_option(self.__instrument_options.layout, 'full')

        if use_cached_bench and layout in SubaruPFI.bench_cache:
            logger.info(f"Getting the bench from cache.")
            bench = SubaruPFI.bench_cache[layout]
        else:
            if layout == 'full':
                raise NotImplementedError('Full bench configuration is no longer supported, use a calibration product.')
            elif layout == 'calibration':
                bench = self.__create_configured_bench()
            else:
                raise NotImplementedError()

            SubaruPFI.bench_cache[layout] = bench

        return bench
    
    def __create_configured_bench(self):
        """
        Create a bench object with real instrument configuration.
        """
        
        layout = self.__get_instrument_option(self.__instrument_options.layout, 'calibration')
        pfs_instdata_path = self.__get_instrument_option(self.__instrument_options.instdata_path, SubaruPFI.DEFAULT_INSTDATA_PATH)
        pfs_black_dots_path = self.__get_instrument_option(self.__instrument_options.blackdots_path, SubaruPFI.DEFAULT_BLACK_DOTS_PATH)
        cobra_coach_dir = self.__get_instrument_option(self.__instrument_options.cobra_coach_dir, SubaruPFI.DEFAULT_COBRA_COACH_DIR)
        cobra_coach_module_version = self.__get_instrument_option(self.__instrument_options.cobra_coach_module_version, None)
        spectrograph_modules = self.__get_instrument_option(self.__instrument_options.spectrograph_modules, [1, 2, 3, 4])
        black_dot_radius_margin = self.__get_instrument_option(self.__instrument_options.black_dot_radius_margin, 1.0)

        # Create the cobra coach temp directory if it does not exist
        if not os.path.isdir(cobra_coach_dir):
            os.makedirs(cobra_coach_dir, exist_ok=True)
            logger.info(f"Created cobra coach temp directory: {cobra_coach_dir}")

        # This is required by cobraCharmer and there seems to be no other way to pass in the
        # calibration product data directory to the library
        os.environ["PFS_INSTDATA_DIR"] = pfs_instdata_path
        
        cobra_coach = CobraCoach("fpga", loadModel=False, trajectoryMode=True, rootDir=cobra_coach_dir)
        cobra_coach.loadModel(version="ALL", moduleVersion=cobra_coach_module_version)
        calib_model = cobra_coach.calibModel

        self.__fix_calib_model(calib_model)

        # Limit spectral modules
        gfm = self.__get_fiber_map()
        cobra_ids_use = np.array([], dtype=np.uint16)
        for sm in spectrograph_modules:
            cobra_ids_use = np.append(cobra_ids_use, gfm.cobrasForSpectrograph(sm))
        logger.info(f"Using {calib_model.nCobras} cobras from {len(spectrograph_modules)} spectrograph modules.")
        logger.info(f"Number of operational cobras: {(calib_model.status == PFIDesign.COBRA_OK_MASK).sum()}.")
        logger.info(f"Number of broken cobras: {((calib_model.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0).sum()}.")
        logger.info(f"Number of broken fibers: {((calib_model.status & PFIDesign.FIBER_BROKEN_MASK) != 0).sum()}.")

        # Set Bad Cobra status for unused spectral modules
        for cobra_id in range(calib_model.nCobras):
            if cobra_id not in cobra_ids_use:
                calib_model.status[cobra_id] = ~PFIDesign.COBRA_OK_MASK
                logger.info(f"Excluding cobraId {cobra_id}.")
        
        # Get the black dots calibration product
        black_dots_calibration_product = BlackDotsCalibrationProduct(pfs_black_dots_path)

        bench = Bench(
            layout=layout,
            calibrationProduct=calib_model,
            # cobraCoach=cobra_coach,
            blackDotsCalibrationProduct=black_dots_calibration_product,
            blackDotsMargin=black_dot_radius_margin,
        )

        return bench, cobra_coach, calib_model
    
    def __fix_calib_model(self, calib_model):
        """
        Implemnents some calib model fixes that are not required anymore.
        """

        ignore_calibration_errors = self.__get_instrument_option(self.__instrument_options.ignore_calibration_errors, False)

        if not ignore_calibration_errors:
            # Run a few quick tests which were required in the past
            zero_centers = calib_model.centers == 0
            zero_link_lengths = (calib_model.L1 == 0) | (calib_model.L2 == 0)
            too_long_link_lengths = (calib_model.L1 > 50) | (calib_model.L2 > 50)

            assert zero_centers.sum() == 0
            assert zero_link_lengths.sum() == 0
            assert too_long_link_lengths.sum() == 0

            # Cobras with broken motor, we don't know where they are
            bad_cobras = (calib_model.status & (calib_model.COBRA_BROKEN_PHI_MASK |
                                          calib_model.COBRA_BROKEN_THETA_MASK |
                                          calib_model.COBRA_BROKEN_MOTOR_MASK)) != 0
            
            assert all(calib_model.phiIn[bad_cobras] == 0.0)
            assert all(calib_model.phiOut[bad_cobras] == -np.pi)
            assert all(calib_model.tht0[bad_cobras] == 0.0)
            assert all(calib_model.tht1[bad_cobras] == 0.0)

            # Cobras with broken fibers, these move but cannot be tracked by the MCS
            # We have a position but we're not sure they're there
            invisible_cobras = (calib_model.status & calib_model.COBRA_INVISIBLE_MASK) != 0

            # TODO: verify?

            broken_fibers = (calib_model.status & calib_model.FIBER_BROKEN_MASK) != 0

            # TODO: verify?

        # bad_cobras = calib_model.status != calib_model.COBRA_OK_MASK

        # # Set some dummy center positions and phi angles for those cobras that have zero centers
        # if bad_cobras.any():
        #     msg = f"Bad cobras: {np.sum(bad_cobras)}"
        #     # calib_model.centers[zero_centers] = np.arange(np.sum(zero_centers)) * 300j
        #     calib_model.phiIn[bad_cobras] = -np.pi
        #     calib_model.phiOut[bad_cobras] = 0
        #     calib_model.tht0[bad_cobras] = 0
        #     calib_model.tht1[bad_cobras] = (2.1 * np.pi) % (2 * np.pi)
        #     logger.warning(msg)

    def get_cobra_centers(self):
        centers = np.array([self.__bench.cobras.centers.real, self.__bench.cobras.centers.imag]).T
        radec, mask = self.__projection.pixel_to_world(centers)

        return centers, radec

    def radec_to_altaz(self,
                       ra, dec, posang, obs_time,
                       epoch=2000.0, pmra=0.0, pmdec=0.0, parallax=1e-6, rv=0.0):
        
        """
        Convert equatorial coordinates to alt-az coordinates at the location of
        the Subaru Telescope.

        Parameters
        ----------
        ra : array of float
            Right Ascension of the targets in degrees
        dec : array of float
            Declination of the targets in degrees
        posang : array of float
            Position angle of the targets in degrees
        obs_time : Time
            Observation time in UTC
        epoch : float
            Epoch of the input coordinates in decimal year. Will be used
            to calculate the position of the targets at ˙obs_time` when
            proper motion is available.
        pmra : array of float
            Proper motion in the RA direction in mas / yr
        pmdec : array of float
            Proper motion in the Dec direction in mas / yr
        parallax : array of float
            Parallax in mas. Due to bugs in the underlying libraries, must
            be non-zero.
        rv : array of float
            Radial velocity in km / s.
        """

        # TODO: move to telescope class

        az, el, inr = DCoeff.radec_to_subaru(ra, dec, posang, obs_time,
                                             epoch=epoch, pmra=pmra, pmdec=pmdec, par=1e-6)

        return az, el, inr
    
    def get_visibility(self, ra, dec, posang, obs_time):
        """
        Check if a target is visible at the Subaru Telescope.
        """

        # TODO: move to telescope class
        
        # Convert pointing center into azimuth and elevation (+ instrument rotator angle)
        alt, az, inr = self.radec_to_altaz(ra, dec, posang, obs_time)

        # Calculate airmass from AZ
        airmass = 1.0 / np.sin(np.radians(90.0 - alt))

        return az > 0, airmass
    
    def __get_reshaped_cobra_config(self, ndim, cobraidx):
        """
        Reshape the cobra configuration to match the shape of the input data.
        """
        
        batch_shape = (Ellipsis,) + (ndim - cobraidx.ndim) * (None,)
        centers = self.__bench.cobras.centers[cobraidx][batch_shape]
        bad_cobra = self.__bench.cobras.hasProblem[cobraidx][batch_shape]
        # bad_cobra = ~self.__bench.cobras.isGood[cobraidx][batch_shape]
        L1 = self.__bench.cobras.L1[cobraidx][batch_shape]
        L2 = self.__bench.cobras.L2[cobraidx][batch_shape]

        return batch_shape, bad_cobra, centers, L1, L2
    
    def __theta_to_local(self, theta, local, cobraidx=None, copy=False):
        """
        Convert theta angles to local, which means theta is measured from the CCW hard stops,
        otherwise theta angles are in global coordinate and phi angles are
        measured from the theta arms. The parameter `local` defines the reference
        system of the input angles.
        """

        cobraidx = cobraidx if cobraidx is not None else ()

        if local:
            return np.copy(theta) if copy else theta
        else:
            return (theta - self.__calib_model.tht0[cobraidx]) % (np.pi * 2)
        
    def __phi_to_local(self, phi, local, cobraidx=None, copy=False):
        cobraidx = cobraidx if cobraidx is not None else ()

        if local:
            return np.copy(phi) if copy else phi
        else:
            return phi - self.__calib_model.phiIn[cobraidx] - np.pi

    def __get_theta_range(self, cobraidx=None):
        cobraidx = cobraidx if cobraidx is not None else ()

        theta_range = (self.__calib_model.tht1 - self.__calib_model.tht0 + np.pi) % (np.pi * 2) + np.pi
        return theta_range[cobraidx]
    
    def __get_phi_range(self, cobraidx=None):
        cobraidx = cobraidx if cobraidx is not None else ()

        phi_range = self.__calib_model.phiOut - self.__calib_model.phiIn
        return phi_range[cobraidx]
    
    def radec_to_fp_pos(self,
                        ra, dec,
                        pointing,
                        epoch=2016.0, pmra=None, pmdec=None, parallax=None, rv=None):
        
        """
        Convert celestial coordinate to focal plane positions.

        Parameters
        ----------
        ra : array of float
            Right Ascension of the targets in degrees
        dec : array of float
            Declination of the targets in degrees
        pointing : Pointing
            Boresight of the telescope and observation time
        epoch : float
            Epoch of the input coordinates in decimal year
        pmra : array of float
            Proper motion in the RA direction in mas / yr
        pmdec : array of float
            Proper motion in the Dec direction in mas / yr
        parallax : array of float
            Parallax in mas
        rv : array of float
            Radial velocity in km / s
        """
        
        cent = np.array([[ pointing.ra ], [ pointing.dec ]])
        pa = pointing.posang

        # The proper motion of the targets used for sky_pfi transformation.
        # The unit is mas/yr, the shape is (2, N)
        if pmra is not None and pmdec is not None:
            pm = np.stack([ pmra, pmdec ], axis=0)
        else:
            pm = None

        # The parallax of the coordinates used for sky_pfi transformation.
        # The unit is mas, the shape is (1, N)
        if parallax is not None:
            parallax = parallax.copy()
        else:
            parallax = None

        # Set pms and parallaxes to 0 if any of them is Nan
        if pm is not None and parallax is not None:
            mask = np.any(np.isnan(pm), axis=0) | np.isnan(parallax)
            pm[:, mask] = 0.0
            parallax[mask] = 1e-7

        # Observation time UTC in format of %Y-%m-%d %H:%M:%S
        obs_time = pointing.obs_time.to_value('iso')
        
        # Input coordinates. Namely. (Ra, Dec) in unit of degree for sky with shape (2, N)
        coords = np.stack([np.atleast_1d(ra), np.atleast_1d(dec)], axis=-1)
        xyin = coords.T

        fp_pos = CoordinateTransform(xyin=xyin,
                                     mode="sky_pfi",
                                     # za=0.0, inr=0.0,     # These are overriden by function
                                     cent=cent,
                                     pa=pa,
                                     pm=pm,
                                     par=parallax,
                                     time=obs_time,
                                     epoch=epoch)
        
        xy = fp_pos[:2, :].T
        return xy[..., 0] + 1j * xy[..., 1]
    
    def fp_pos_to_rel_fp_pos(self, fp_pos, cobraidx):
        """
        Convert focal plane positions to relative focal plane positions. It also returns
        a few precomputed parameters.

        Parameters
        ----------
        fp_pos : np.ndarray
            Focal plane positions with shape (N,) where N is the number of targets.
        cobraidx : np.ndarray
            Cobra indices for each focal plane position with shape (N,) where N is the number of targets.

        Returns
        -------
        batch_shape : tuple
            The shape of the batch.
        bad_cobra : np.ndarray
            Boolean array indicating which cobras are in a bad configuration.
        centers : np.ndarray
            The centers of the cobras.
        L1 : np.ndarray
            The lengths of the first arm of the cobras.
        L2 : np.ndarray
            The lengths of the second arm of the cobras.
        d : np.ndarray
            Distance from the main axis of the cobra
        d_2 : np.ndarray
            Square of the distance from the main axis of the cobra
        """

        batch_shape, bad_cobra, centers, L1, L2 = self.__get_reshaped_cobra_config(fp_pos.ndim, cobraidx)

        rel_fp_pos = fp_pos - centers
        d = np.abs(rel_fp_pos)
        d_2 = d * d

        return batch_shape, bad_cobra, rel_fp_pos, centers, L1, L2, d, d_2 
        
    def fp_pos_to_cobra_angles(self, fp_pos, cobraidx):
        """
        Convert focal plane positions to local cobra angles, which are measured from
        the cobra hard stops.

        For each focal plane positions, two sets of theta and phi angles are calculated

        The function is a slightly improved version of the cobraCharmer.PFI.positionsToAngles
        to allow calculating the angles for multiple targets at the same time.

        Parameters
        ----------
        fp_pos : np.ndarray
            Focal plane positions with shape (N,) where N is the number of targets.
        cobraidx : np.ndarray
            Cobra indices for each focal plane position with shape (N,) where N is the number of targets.

        Returns
        -------
        theta : np.ndarray
            Theta angles with shape (N, 2) where N is the number of targets. The first column is always filled,
            the second column is filled only if secondary solutions are available.
        phi : np.ndarray
            Phi angles with shape (N, 2) where N is the number of targets. The first column is always filled,
            the second column is filled only if secondary solutions are available.
        d : np.ndarray
            Distance from the main axis of the cobra
        eb_pos : np.ndarray
            Elbow positions with shape (N, 2) where N is the number of targets.
        flags : np.ndarray
            Flags with shape (N, 2) where N is the number of targets indicating the validity of the solutions.
        """

        # TODO: this could be taken from cobraOps but this is a better vectorized version of
        #       what's available there

        batch_shape, bad_cobra, rel_fp_pos, centers, L1, L2, d, d_2 = \
            self.fp_pos_to_rel_fp_pos(fp_pos, cobraidx)

        L1_2 = L1 ** 2
        L2_2 = L2 ** 2

        arp = np.angle(rel_fp_pos)

        # This might generate warnings because we don't filter on unreachable targets yet
        with np.errstate(invalid='ignore'):
            ang1 = np.arccos((L1_2 + L2_2 - d_2) / (2 * L1 * L2))
            ang2 = np.arccos((L1_2 + d_2 - L2_2) / (2 * L1 * d))
        
        phi_in = self.__bench.cobras.phiIn[cobraidx][batch_shape] + np.pi
        phi_out = self.__bench.cobras.phiOut[cobraidx][batch_shape] + np.pi
        theta0 = self.__bench.cobras.tht0[cobraidx][batch_shape]
        theta1 = self.__bench.cobras.tht1[cobraidx[batch_shape]]
        
        # Return values
        phi = np.full(fp_pos.shape + (2,), np.nan)
        theta = np.full(fp_pos.shape + (2,), np.nan)
        flags = np.full(fp_pos.shape + (2,), 0, dtype=int)

        # Handle the various cases of problematic configurations
        
        # - too far away, return theta= spot angle and phi=PI
        too_far_from_center = ~bad_cobra & (d > L1 + L2)
        flags[too_far_from_center, 0] |= CobraAngleFlags.TOO_FAR_FROM_CENTER
        phi[too_far_from_center, 0] = np.pi
        theta[too_far_from_center, 0] = (arp - theta0)[too_far_from_center] % (2 * np.pi)
        
        # - too close to center, theta is undetermined, return theta=spot angle and phi=0
        too_close_to_center = ~bad_cobra & (d < np.abs(L1 - L2))
        flags[too_close_to_center, 0] |= CobraAngleFlags.TOO_CLOSE_TO_CENTER
        phi[too_close_to_center, 0] = 0.0
        theta[too_close_to_center, 0] = (arp - theta0)[too_close_to_center] % (2 * np.pi)

        # -- Solution type 1
        
        # - the regular solutions, phi angle is between 0 and pi, no checking for phi hard stops
        solution_ok = ~bad_cobra & ~too_far_from_center & ~too_close_to_center
        flags[solution_ok, 0] |= CobraAngleFlags.SOLUTION_OK
        phi[solution_ok, 0] = (ang1 - phi_in)[solution_ok]
        theta[solution_ok, 0] = (arp + ang2 - theta0)[solution_ok] % (2 * np.pi)
        # theta in overlap region
        in_overlapping_region = solution_ok & (theta[..., 0] <= (theta1 - theta0) % (2 * np.pi))
        flags[in_overlapping_region, 0] |= CobraAngleFlags.IN_OVERLAPPING_REGION
        
        # -- Solution type 2

        solution_ok_2 = solution_ok & (ang1 <= np.pi / 2) & (ang1 > 0)
        phi_ok_2 = (phi_in <= -ang1)
        flags[solution_ok_2 & phi_ok_2, 1] |= CobraAngleFlags.SOLUTION_OK
        flags[solution_ok_2, 1] |= CobraAngleFlags.PHI_NEGATIVE
        # phiIn < 0
        phi[solution_ok_2, 1] = (-ang1 - phi_in)[solution_ok_2]
        theta[solution_ok_2, 1] = (arp - ang2 - theta0)[solution_ok_2] % (2 * np.pi)
        # theta in overlap region
        in_overlapping_region = solution_ok_2 & (theta[..., 1] <= (theta1 - theta0) % (2 * np.pi))
        flags[in_overlapping_region, 1] |= CobraAngleFlags.IN_OVERLAPPING_REGION

        # -- Solution type  3

        solution_ok_3 = solution_ok & (ang1 > np.pi / 2) & (ang1 < np.pi)
        phi_ok_3 = (phi_out >= 2 * np.pi - ang1)
        flags[solution_ok_3 & phi_ok_3, 1] |= CobraAngleFlags.SOLUTION_OK
        flags[solution_ok_3, 1] |= CobraAngleFlags.PHI_BEYOND_PI
        # phiOut > np.pi
        phi[solution_ok_3, 1] = (2 * np.pi - ang1 - phi_in)[solution_ok_3]
        theta[solution_ok_3, 1] = (arp - ang2 - theta0)[solution_ok_3] % (2 * np.pi)
        # theta in overlap region
        in_overlapping_region = solution_ok_3 & (theta[..., 1] <= (theta1 - theta0) % (2 * np.pi))
        flags[in_overlapping_region, 1] |= CobraAngleFlags.IN_OVERLAPPING_REGION

        # Check out of range angles

        theta_range = (theta1 - theta0 + np.pi) % (2 * np.pi) + np.pi
        theta_out_of_range = solution_ok & ((theta[..., 0] < 0) | (theta_range < theta[..., 0]))
        flags[theta_out_of_range, 0] |= CobraAngleFlags.THETA_OUT_OF_RANGE
        theta_out_of_range = (solution_ok_2 | solution_ok_3) & ((theta[..., 1] < 0) | (theta_range < theta[..., 1]))
        flags[theta_out_of_range, 1] |= CobraAngleFlags.THETA_OUT_OF_RANGE

        phi_range = phi_out - phi_in
        phi_out_of_range = solution_ok & ((phi[..., 0] < 0) | (phi_range < phi[..., 0]))
        flags[phi_out_of_range, 0] |= CobraAngleFlags.PHI_OUT_OF_RANGE
        phi_out_of_range = (solution_ok_2 | solution_ok_3) & ((phi[..., 1] < 0) | (phi_range < phi[..., 1]))
        flags[phi_out_of_range, 1] |= CobraAngleFlags.PHI_OUT_OF_RANGE

        # Calculate the elbow positions while we're at it
        eb_pos = np.full_like(theta, np.nan, dtype=complex)
        eb_pos[solution_ok, 0] = (centers + L1 * np.exp(1j * (theta[..., 0] + theta0)))[solution_ok]
        eb_pos[solution_ok_2 | solution_ok_3, 1] = (centers + L1 * np.exp(1j * (theta[..., 1] + theta0)))[solution_ok_2 | solution_ok_3]
                    
        return theta, phi, d, eb_pos, flags
    
    def cobra_angles_to_fp_pos(self, theta, phi, cobraidx):
        """
        Convert local cobra angles to focal plane position. The angles are measured
        from the cobra hard stops.
        """

        batch_shape, bad_cobra, centers, L1, L2 = self.__get_reshaped_cobra_config(theta.ndim, cobraidx)

        phi_in = self.__bench.cobras.phiIn[cobraidx][batch_shape]
        phi_out = self.__bench.cobras.phiOut[cobraidx][batch_shape]
        theta0 = self.__bench.cobras.tht0[cobraidx][batch_shape]
        theta1 = self.__bench.cobras.tht1[cobraidx[batch_shape]]

        # Do some sanity checks
        # TODO: move this to Netflow instead

        # TODO: use __get_theta_range instead?
        theta_range = (theta1 - theta0 + np.pi) % (2 * np.pi) + np.pi
        if np.any(0 > theta) or np.any(theta_range < theta):
            logger.error('Some theta angles are out of range')
        
        phi_range = phi_out - phi_in
        if np.any(0 > phi) or np.any(phi_range < phi):

            # TODO even when the phi angles are set to home position, this
            #      is triggered for cobra 46 an 49
            #      in both cases phi_in = 0 and phi_out = -pi
            #      these cobras are marked as broken

            logger.error('Some phi angles are out of range')

        ang1 = theta + theta0
        ang2 = ang1 + phi + phi_in
        fp_pos = centers + L1 * np.exp(1j * ang1) + L2 * np.exp(1j * ang2)

        return fp_pos
    
    def create_cobra_state(self, cobraidx):
        """
        Create a new cobra state object.
        """

        cobra_state = SimpleNamespace()

        cobra_state.mode = 'normal'                 # choices: 'normal', 'theta', 'phi', only 'normal' is implemented
        cobra_state.constant_speed = False          # self.__cobra_coach.constantSpeedMode
        cobra_state.max_segments = self.__cobra_coach.maxSegments
        cobra_state.max_steps_per_segment = self.__cobra_coach.maxStepsPerSeg
        cobra_state.max_total_steps =  2000         # self.__cobra_coach.maxTotalSteps
        cobra_state.two_step_method = True          # when True, make half-steps
        cobra_state.phi_threshold = 0.6             # defined in cobra_coach.calculus.phiThreshold, minimum phi angle to move in the beginning
        cobra_state.time_step = 20                  # defined in trajectories, the trajectory time resolution in step units

        cobra_state.tolerance = 0.005               # tolerance for target positions in pixels
        cobra_state.fast_move_threshold = 20.0      # the threshold for using slow or fast motor maps
        
        cobra_state.theta_margin = np.deg2rad(15.0),  # the minimum theta angles to the theta hard stops

        # Cobras being simulated. Any number of targets per cobra is supported.
        cobra_state.cobraidx = cobraidx

        # Cobras with broken motors
        cobra_state.broken = ((self.__calib_model.status[cobraidx] & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0)

        cobra_state.theta_range = self.__get_theta_range(cobraidx)
        cobra_state.phi_range = self.__get_phi_range(cobraidx)

        # Current positions
        cobra_state.theta = np.full_like(cobraidx, np.nan, dtype=float)
        cobra_state.phi = np.full_like(cobraidx, np.nan, dtype=float)
        cobra_state.fp_pos = np.full_like(cobraidx, np.nan + 1j * np.nan, dtype=complex)

        return cobra_state
    
    def create_trajectories(self, cobra_state):
        """
        Create a new trajectories object.
        """

        trajectories = SimpleNamespace()

        # Cobras being simulated. Any number of targets per cobra is supported.
        trajectories.cobraidx = cobra_state.cobraidx

        trajectories.time_step = cobra_state.time_step

        trajectories.steps = 0
        trajectories.mask = []
        trajectories.theta = []
        trajectories.phi = []
        trajectories.fp_pos = []                # Tip focal plane position
        trajectories.eb_pos = []                # Elbow focal plane position

        return trajectories
    
    def create_moves(self, cobra_state, tries):
        """
        Create a new object to keep track of cobra movements.
        """

        moves_shape = cobra_state.cobraidx.shape + (tries,)
        moves = SimpleNamespace()

        moves.fp_pos = []
        moves.theta = []
        moves.theta_steps = []
        moves.theta_ontime = []
        moves.theta_fast = []
        moves.phi = []
        moves.phi_steps = []
        moves.phi_ontime = []
        moves.phi_fast = []

        return moves
    
    def set_cobra_state_to_home(self,
                                cobra_state, /,
                                theta_enabled=True,
                                phi_enabled=True,
                                theta_ccw=True,):
                                
        """
        Move the cobras to the home position.

        Parameters
        ----------
        theta_enabled : bool
            Enable theta motor or not.
        theta_ccw : bool
            Move theta motor in CCW direction when True.
        phi_enabled : bool
            Enable phi motor or not.
        """

        # If not specified, create a new cobra state with cobras reset to home
        # position. This state is stored within the cobraInfo variable of cobra coach
        # but we simplified it here for the sake of simulation.

        # Calculate the starting angles. These are local angles measured from the CCW hard stops.
        if theta_enabled:
            if theta_ccw:
                cobra_state.theta = np.zeros_like(cobra_state.cobraidx, dtype=float)
            else:
                cobra_state.theta = self.__get_theta_range(cobra_state.cobraidx)

        if phi_enabled:
            cobra_state.phi = np.zeros_like(cobra_state.cobraidx, dtype=float)

        # Calculate the new focal plane positions
        cobra_state.fp_pos = self.cobra_angles_to_fp_pos(cobra_state.theta, cobra_state.phi, cobra_state.cobraidx)
    
    def simulate_trajectories(self,
                              cobra_state,
                              theta_target, phi_target,
                              use_scaling=True,
                              max_segments=1,
                              max_total_steps=2000):
        """
        Simulate cobra movement trajectories.
        """

        # This function is roughly the equivalent of engineer.createTrajectory

        # TODO: bring out as parameters, or even better, add them to the cobra_state object
        tries = 8

        # Create an object to keep track of the trajectories
        trajectories = self.create_trajectories(cobra_state)

        # Create object to keep track of the movements
        moves = self.create_moves(cobra_state, tries)

        if cobra_state.two_step_method:
            # Perform movement in two steps

            # Limit phi angle for first two tries
            phi_limit = np.pi / 3 - self.__calib_model.phiIn[cobra_state.cobraidx] - np.pi
            theta_via = np.copy(theta_target)
            phi_via = np.copy(phi_target)

            # Constrain movement within the range of the motors
            phi_out_of_range = phi_target > phi_limit
            phi_via = np.where(phi_out_of_range, phi_limit, phi_via)
            theta_via = np.where(phi_out_of_range, theta_target + (phi_target - phi_limit) / 2, theta_via)
            
            theta_range = self.__get_theta_range(cobra_state.cobraidx)
            theta_out_of_range = phi_out_of_range & (theta_via > theta_range)
            theta_via = np.where(theta_out_of_range, theta_range, theta_via)

            # Make the first step

            self.__move_theta_phi(cobra_state,
                                  trajectories,
                                  moves,
                                  theta_via, phi_via,
                                  relative=False,
                                  local=True,
                                  tries=2,
                                  theta_fast=True,                  # theta moves fast in first step
                                  phi_fast=True,
                                  use_scaling=False,
                                  max_segments=2 * max_segments,
                                  max_total_steps=2 * max_total_steps)

            # Make the second step

            self.__move_theta_phi(cobra_state,
                                  trajectories,
                                  moves,
                                  theta_target, phi_target,
                                  relative=False,
                                  local=True,
                                  tries=2,
                                  theta_fast=False,                 # theta moves slowly in second step
                                  phi_fast=True,
                                  use_scaling=use_scaling,
                                  max_segments=max_segments,
                                  max_total_steps=max_total_steps)
        else:
            # Perform movements in a single step

            self.__move_theta_phi(cobra_state,
                                  trajectories,
                                  moves,
                                  theta_target, phi_target,
                                  relative=False,
                                  local=True,
                                  tries=tries,
                                  theta_fast=False,                 # theta moves slowly
                                  phi_fast=True,
                                  use_scaling=use_scaling,
                                  max_segments=max_segments,
                                  max_total_steps=max_total_steps)

    def __move_theta_phi(self,
                         cobra_state,
                         trajectories,
                         moves,
                         theta, phi,
                         relative=False,    # angles are offsets from current positions or not
                         local=True,        # the angles are from the CCW hard stops or not
                         tries=6,
                         theta_fast=False, phi_fast=False, # use fast motor maps
                         use_scaling=True, 
                         max_segments=1, 
                         max_total_steps=1
                         ):
        
        """
        Perform a full move from the current positions, as defined in cobra_state, to the
        given theta and phi angles. The angles can be relative to the current positions, in case
        `relative` is True, or they can be measured from the CCW hard stops, in case `local` is True.

        If local is True, both theta and phi angles are measured from CCW hard stops,
        otherwise theta angles are in global coordinate and phi angles are
        measured from the phi arms.

        This function is the equivalent of engineer.moveThetaPhi

        Parameters
        ----------
        cobra_state : SimpleNamespace
            The cobra state object.
        trajectories : SimpleNamespace
            Object to store the trajectory segments.
        moves : SimpleNamespace
            Objects to store information about the cobra movements.
        theta : np.ndarray
            The target theta angles.
        phi : np.ndarray
            The target phi angles.
        relative : bool
            True if the angles are relative to current positions.
        local : bool
            True if the angles are measured from the CCW hard stops,
            otherwise theta angles are in global coordinate and phi angles are
            measured from the theta arms.
        tries : int
            The number of tries to converge.
        theta_fast : bool
            Use fast motor maps for theta. Threshold is set in `cobra_state`.
        phi_fast : bool
            Use fast motor maps for phi. Threshold is set in `cobra_state`.
        """
        
        theta_range = self.__get_theta_range(cobra_state.cobraidx)
        
        theta_fast = np.full(cobra_state.cobraidx.shape, theta_fast, dtype=bool)
        phi_fast = np.full(cobra_state.cobraidx.shape, phi_fast, dtype=bool)

        shape = cobra_state.cobraidx.shape
        done_mask = np.full(shape, False, dtype=bool)           # Cobra has converged
        not_done_mask = np.full(shape, True, dtype=bool)        # Cobra hasn't converged yet
        far_away_mask = np.full(shape, True, dtype=bool)        # Cobra is far away from the target position
        theta_at = np.full(shape, False, dtype=float)           # Current theta angles
        phi_at = np.full(shape, False, dtype=float)             # Current phi angles

        if relative:
            # Target angles are relative to current positions
            theta_target = cobra_state.theta + theta
            phi_target = cobra_state.phi + phi
        else:
            # Target angles are absolute, convert them to local
            theta_target = self.__theta_to_local(theta, local, cobraidx=cobra_state.cobraidx, copy=True)
            phi_target = self.__phi_to_local(phi, local, cobraidx=cobra_state.cobraidx, copy=True)

            if not local:
                # Theta is outside the margins
                mask = theta_target < cobra_state.theta_margin
                theta_target[mask] += np.pi * 2

                # Theta is out of range
                mask = theta_target > theta_range - cobra_state.theta_margin
                theta_target[mask] = theta_range[mask] - cobra_state.theta_margin

                # TODO: shouldn't this be reported as a warning?
                #       this cobra will never converge

        # Calculate the target focal plane positions
        fp_pos_target = self.cobra_angles_to_fp_pos(theta_target, phi_target, cobra_state.cobraidx)

        # In simulation mode, cc.camResetStack is a no-op
        # In the original code, cobras are sent to home position, here we assume that
        # cobra_state is already set to home position, or any position that makes sense

        for j in range(tries):
            # If the cobras are far away enough from the target position, we can move them fast
            theta_fast[not_done_mask] = far_away_mask[not_done_mask] & theta_fast[not_done_mask]
            phi_fast[not_done_mask] = far_away_mask[not_done_mask] & phi_fast[not_done_mask]

            theta_scaling = ~far_away_mask
            phi_scaling = ~far_away_mask

            # Perform the movement. This will update cobra_state to reflect the new positions
            theta_at[not_done_mask], phi_at[not_done_mask] = self.__move_to_angles(cobra_state,
                                                                                   trajectories,
                                                                                   moves,
                                                                                   theta_target, phi_target,
                                                                                   theta_fast, phi_fast,
                                                                                   theta_scaling=theta_scaling,
                                                                                   phi_scaling=phi_scaling,
                                                                                   local=True,
                                                                                   mask=not_done_mask)
            
            # Keep track of the positions
            fp_pos_at = cobra_state.fp_pos
            d = np.abs(fp_pos_at - fp_pos_target)
            done_mask[d < cobra_state.tolerance] = True
            newly_done_mask = done_mask & not_done_mask
            far_away_mask = (d > cobra_state.fast_move_threshold) & far_away_mask
            
            # Keep track of the positions
            # TODO: some of these should be returned from __move_to_angles
            #       instead of the positions because those are available in cobra_state
            if moves is not None:
                moves.fp_pos.append(np.copy(cobra_state.fp_pos))
                moves.theta.append(np.copy(cobra_state.theta))
                # moves.theta_steps = []
                # moves.theta_ontime = []
                moves.theta_fast.append(np.copy(theta_fast))
                moves.phi.append(np.copy(cobra_state.phi))
                # moves.phi_steps = []
                # moves.phi_ontime = []
                moves.phi_fast.append(np.copy(phi_fast))

            # Keep track cobras which haven't converged yet
            if np.any(newly_done_mask):
                not_done_mask &= ~newly_done_mask

            # Check if all cobras have converged
            if not np.any(not_done_mask):
                break

        if np.any(not_done_mask):
            logger.warning(f"Failed to converge for {np.sum(not_done_mask)} cobras.")

        return moves

    def __move_to_angles(self,
                         cobra_state,
                         trajectories,
                         moves,
                         theta_target, phi_target,
                         theta_fast, phi_fast,
                         theta_scaling=True, phi_scaling=True,
                         local=True,
                         mask=None):
        
        """
        Move cobras to the given theta and phi angles.

        If local is True, both theta and phi angles are measured from CCW hard stops,
        otherwise theta angles are in global coordinate and phi angles are
        measured from the phi arms.

        This function is closely based on CobraCoach.moveToAngles but many of the checks
        are removed for simplicity.
        """

        if theta_target is not None:
            theta_to = self.__theta_to_local(theta_target, local, cobraidx=cobra_state.cobraidx)
            theta_from = cobra_state.theta
            theta_delta = theta_to - theta_from
        else:
            theta_to = None
            theta_from = None
            theta_delta = None

        if phi_target is not None:
            phi_to = self.__phi_to_local(phi_target, local, cobraidx=cobra_state.cobraidx)
            phi_from = cobra_state.phi
            phi_delta = phi_to - phi_from
        else:
            phi_to = None
            phi_from = None
            phi_delta = None

        # Calculate the trajectory along the movement
        self.__move_delta_angles(cobra_state,
                                 trajectories,
                                 moves,
                                 theta_delta, phi_delta,
                                 theta_fast, phi_fast,
                                 theta_scaling=theta_scaling,
                                 phi_scaling=phi_scaling,
                                 mask=mask)

    def __move_delta_angles(self,
                            cobra_state,
                            trajectories,
                            moves,
                            theta_delta, phi_delta,
                            theta_fast, phi_fast,
                            theta_scaling, phi_scaling,
                            mask=None):
        """
        Move cobras by the given amount.

        This is a simplified version of CobraCoach.moveDeltaAngles with most of the checks removed.
        """

        # We do not implement theta-only and phi-only modes here, only normal mode
        
        # Current cobra angles
        theta_from = cobra_state.theta
        phi_from = cobra_state.phi

        if not cobra_state.constant_speed:
            theta_steps, phi_steps, theta, phi = \
                self.__calculate_steps(cobra_state,
                                       theta_delta, phi_delta,
                                       theta_from, phi_from,
                                       theta_fast, phi_fast,
                                       max_steps = cobra_state.max_total_steps)

            self.__move_steps(cobra_state,
                              trajectories,
                              moves,
                              theta_steps, phi_steps, 
                              theta_expected=theta, phi_expected=phi,
                              theta_fast=theta_fast, phi_fast=phi_fast,
                              mask=mask)
        else:
            raise NotImplementedError()
        
            # # Determine the maximum number of required segments
            # segments = 0

            # theta_mm = np.where(theta_fast, self.__cobra_coach.mmTheta[cobra_state.cobraidx], self.__cobra_coach.mmThetaSlow[cobra_state.cobraidx])
            # phi_mm = np.where(phi_fast, self.__cobra_coach.mmPhi[cobra_state.cobraidx], self.__cobra_coach.mmPhiSlow[cobra_state.cobraidx])

            # theta_n = self.__calculate_segments(cobra_state, theta_delta, theta_from, theta_mm)
            # phi_n = self.__calculate_segments(cobra_state, phi_delta, phi_from, phi_mm)

            # n = np.max([np.max(theta_n), np.max(phi_n)])
            # n = np.min([n, cobra_state.max_segments])

    def __move_steps(self,
                     cobra_state,
                     trajectories,
                     moves,
                     theta_steps, phi_steps,
                     theta_fast=False, phi_fast=False,
                     theta_expected=None, phi_expected=None,
                     theta_ontime=None, phi_ontime=None,
                     theta_speed=None, phi_speed=None,
                     force=False,
                     n_segments=0,
                     mask=None):
        """
        
        This function is based on cobraCoach.moveSteps but most of the check are removed and
        only the trajectory simulation branch is implemented.

        Parameters
        ----------
        cobra_state : SimpleNamespace
            Current cobra state.
        trajectories : SimpleNamespace
            Object to store the trajectory segments.
        theta_steps : np.ndarray of int
            Number of steps for theta motor.
        phi_steps : np.ndarray of int
            Number of steps for phi motor.
        theta_fast : np.ndarray of bool
            Use fast theta motor maps.
        phi_fast : np.ndarray of bool
            Use fast phi motor maps.
        theta_expected : np.ndarray of float
            Expected theta angles.
        phi_expected : np.ndarray of float
            Expected phi angles.
        theta_ontime : np.ndarray of float
            On-time theta angles.
        phi_ontime : np.ndarray of float
            On-time phi angles.
        theta_speed : np.ndarray of float
            Speed of theta motor.
        phi_speed : np.ndarray of float
            Speed of phi motor.
        force : bool
            Force the move.
        n_segments : int
            Number of segments.
        mask : np.ndarray of bool
            Mask of cobras to move.
        """

        # TODO: we don't need this, everything should be limited by max_total_steps
        # theta_steps = np.clip(theta_steps, -10000, 10000)
        # phi_steps = np.clip(phi_steps, -6000, 6000)

        if n_segments == 0:
            # Interpolate the move into finite increments
            theta_tr, phi_tr = self.__interpolate_moves(cobra_state,
                                                        theta_steps, phi_steps,
                                                        theta_fast, phi_fast,
                                                        mask=mask)
        elif n_segments > 0:
            raise NotImplementedError()
            theta_tr, phi_tr = self.__interpolate_segments(cobra_state,
                                                         theta_steps, phi_steps,
                                                         theta_fast, phi_fast,
                                                         mask=mask)

            # Here comes pfi.moveSteps which is different from this function
            # record the move to the trajectory history
        else:
            raise NotImplementedError()
        
        # Keep track of the trajectories
        if trajectories is not None:
            trajectories.steps += theta_tr.shape[-1]
            trajectories.mask.append(mask)
            trajectories.theta.append(theta_tr)
            trajectories.phi.append(phi_tr)
            trajectories.fp_pos = []
            trajectories.eb_pos = []

        if n_segments == 0:
            self.__pfi_move_steps(cobra_state,
                                  theta_steps, phi_steps,
                                  theta_fast, phi_fast,
                                  theta_ontime, phi_ontime,
                                  mask=mask)
        else:
            for n in range(n_segments):
                self.__pfi_move_steps(cobra_state,
                                      theta_steps[n], phi_steps[n],
                                      theta_ontime=theta_ontime[n],
                                      phi_ontime=phi_ontime[n],
                                      mask=mask)

        # TODO: implement the two other modes?        
        # if theta_mode
        #     pass
        # elif phi mone:
        #     pass
        # else:
        theta = np.copy(cobra_state.theta)
        phi = np.copy(cobra_state.phi)

        # TODO: what are these expected thetas?
        #       they're set to nan if the input is None in the original code
        theta[mask] += theta_expected[mask]
        phi[mask] += phi_expected[mask]

        targets = self.cobra_angles_to_fp_pos(theta, phi, np.arange(theta.size))
        
        
                
        # TODO: continue from cobraCoach.py#L620

        # theta = np.copy(cobra_state.theta)
        # phi = np.copy(cobra_state.phi)
        # theta += theta_expected
        # phi += phi_expected
        # targets = self.pfi.anglesToPositions(self.allCobras, theta, phi)

        # pos[self.visibleIdx] = targets[self.visibleIdx]

        # # update status
        # # ...

    def __pfi_move_steps(self, 
                         cobra_state,
                         theta_steps, phi_steps,
                         theta_steps_wait=None, phi_steps_wait=None,
                         theta_fast=True, phi_fast=True,
                         theta_ontime=None, phi_ontime=None,
                         interval=2.5,
                         force=False,
                         mask=None):
        
        """
        Move cobras with theta and phi steps
        """
        
        # This should be cobraCharmer.PFI.moveSteps

        # No idea if is has any side effect

        # TODO: find where the cobra positions are updated because we might not even need this function at all
        # raise NotImplementedError()

        pass

    def __interpolate_moves(self,
                            cobra_state,
                            theta_steps,
                            phi_steps,
                            theta_fast,
                            phi_fast,
                            mask=None):
        """
        Interpolate a move into step increments. The quantum of increments is `time_step` and
        the increments are distributed over a total number of steps. Then integrate these increments
        and interpolate the motor model to get the angles, including the starting point.

        The starting point of the move is the current cobra positions.

        Parameters
        ----------
        cobra_state : SimpleNamespace
            Current cobra state.
        theta_steps : np.ndarray of int
            Number of steps for theta motor.
        phi_steps : np.ndarray of int
            Number of steps for phi motor.
        theta_fast : np.ndarray of bool
            Use fast theta motor maps.
        phi_fast : np.ndarray of bool
            Use fast phi motor maps.

        Returns
        -------
        theta_tr : np.ndarray of float
            Trajectory of theta motor.
        phi_tr : np.ndarray of float
            Trajectory of phi motor.
        """

        # I guess we start from cobra_state.phi and cobra_state.theta and make a step of
        # theta_steps and phi_steps here, we'll see

        time_step = cobra_state.time_step
        max_steps = np.max(np.abs([theta_steps, phi_steps]))
        size = int(np.ceil(max_steps / time_step))

        # Select the correct theta model for each (virtual) cobra
        theta_model_step, theta_model_offset = self.__get_theta_model(cobra_state, theta_steps, theta_fast)
        phi_model_step, phi_model_offset = self.__get_phi_model(cobra_state, phi_steps, phi_fast)

        # For all (virtual) cobras, calculate the theta and phi step increments
        theta_inc, phi_inc = self.__interpolate_steps(cobra_state, theta_steps, phi_steps, max_steps)

        # Calculate trajectory
        
        # Create the trajectory arrays for each (virtual) cobra and
        # copy the current cobra positions to the 0th item
        shape = cobra_state.cobraidx.shape + (size + 1,)
        
        theta_tr = np.zeros(shape)
        theta_tr[..., 0] = cobra_state.theta
        
        phi_tr = np.zeros(shape)
        phi_tr[..., 0] = cobra_state.phi

        # Interpolate the motor model to the start angle
        theta_start = self.__batch_interp(cobra_state.theta, theta_model_offset, theta_model_step)
        phi_start = self.__batch_interp(cobra_state.phi, phi_model_offset, phi_model_step)
        
        # Integrate the step increments and interpolate the motor model to get the angles
        theta_tr[..., 1:] = self.__batch_interp(np.cumsum(theta_inc, axis=-1) + theta_start[..., None], theta_model_step, theta_model_offset)
        phi_tr[..., 1:] = self.__batch_interp(np.cumsum(phi_inc, axis=-1) + phi_start[..., None], phi_model_step, phi_model_offset)

        return theta_tr, phi_tr

    def __interpolate_steps(self, cobra_state, theta_step, phi_step, max_step):
        """
        Interpolate the steps for the theta and phi motors by taking the step delays into account.
        The steps are made in increments of `cobra_state.time_step` and the increment is calculated
        for each point in time. Some rotations are delayed.

        This function is based on calculus.interpolateSteps but vectorized.

        Parameters
        ----------
        cobra_state : SimpleNamespace
            Current cobra state.
        theta_step : np.ndarray of int
            Total number of steps for theta motor.
        phi_step : np.ndarray of int
            Total number of steps for phi motor.
        max_step : int
            Maximum number of total steps for both motors.

        Returns
        -------
        theta_inc : np.ndarray of int
            Number of steps for theta motor as a function of time.
        phi_inc : np.ndarray of int
            Number of steps for phi motor as a function of time.
        """

        time_step = cobra_state.time_step
        size = int(np.ceil(max_step / time_step))
        batch_shape = cobra_state.cobraidx.shape

        # Depending on the delays, calculate the theta and phi steps
        # This is the same as the original code but vectorized
        
        # time_step is the trajectory time step resolution in step units
        time = np.arange(size) * time_step

        # Set delay parameters for safer operation, same logic as in pfi.py
        theta_delay, phi_delay = self.__get_step_delay(cobra_state, theta_step, phi_step, max_step)

        # Calculate all theta increments

        # shapes:
        # - theta_mask: (batch_shape)
        # - theta_steps: (batch_shape, size)
        theta_sign = np.where(theta_step < 0, -1, 1)

        theta_mask = (time < theta_delay[..., None] + theta_step[..., None] * theta_sign[..., None]) & \
                     (time > theta_delay[..., None] - time_step)
        
        theta_inc = np.zeros(batch_shape + (size,), int)
        theta_inc[theta_mask] = time_step
        
        exceed = (theta_delay[..., None] - time)[theta_mask]
        exceed_mask = exceed > 0
        theta_inc[theta_mask][exceed_mask] = theta_inc[theta_mask][exceed_mask] - exceed[exceed_mask]

        exceed = (time + time_step - theta_delay[..., None] - theta_step[..., None] * theta_sign[..., None])[theta_mask]
        exceed_mask = exceed > 0
        theta_inc[theta_mask][exceed_mask] = theta_inc[theta_mask][exceed_mask] - exceed[exceed_mask]

        assert np.all(np.abs(theta_inc.sum(axis=-1) - theta_step * theta_sign) <= time_step)
        
        # Calculate all phi increments

        phi_sign = np.where(phi_step < 0, -1, 1)

        phi_mask = (time < phi_delay[..., None] + phi_step[..., None] * phi_sign[..., None]) & \
                   (time > phi_delay[..., None] - time_step)
        
        phi_inc = np.zeros(batch_shape + (size,), int)
        phi_inc[phi_mask] = time_step

        exceed = (phi_delay[..., None] - time)[phi_mask]
        exceed_mask = exceed > 0
        phi_inc[phi_mask][exceed_mask] = phi_inc[phi_mask][exceed_mask] - exceed[exceed_mask]

        exceed = (time + time_step - phi_delay[..., None] - phi_step[..., None] * phi_sign[..., None])[phi_mask]
        exceed_mask = exceed > 0
        phi_inc[phi_mask][exceed_mask] = phi_inc[phi_mask][exceed_mask] - exceed[exceed_mask]

        assert np.all(np.abs(phi_inc.sum(axis=-1) - phi_step * phi_sign) <= time_step)

        ####

        # Switch back the sign
        theta_inc *= theta_sign[..., None]
        phi_inc *= phi_sign[..., None]

        return theta_inc, phi_inc

    def __get_theta_model(self, cobra_state, theta_delta, theta_fast):
        """
        Return the calibration model of the theta motor for each (virtual) cobra
        depending on the direction and speed.
        """

        # This is used at at least two locations in the original code:
        # - CobraCoach.calculateSteps
        # - calculus.interpolateMoves

        calib_steps = self.__calib_model.thtOffsets.shape[1]
        batch_shape = np.shape(cobra_state.cobraidx)
        shape = batch_shape + (calib_steps,)

        theta_offset = self.__calib_model.thtOffsets[cobra_state.cobraidx]
        theta_steps = np.full(shape, np.nan)

        # Depending on the direction of the movement and speed, select the correct model

        mask = (theta_delta >= 0) & theta_fast
        theta_steps[mask] = self.__calib_model.posThtSteps[cobra_state.cobraidx[mask], :]

        mask = (theta_delta >= 0) & ~theta_fast
        theta_steps[mask] = self.__calib_model.posThtSlowSteps[cobra_state.cobraidx[mask], :]

        mask = (theta_delta < 0) & theta_fast
        theta_steps[mask] = self.__calib_model.negThtSteps[cobra_state.cobraidx[mask], :]

        mask = (theta_delta < 0) & ~theta_fast
        theta_steps[mask] = self.__calib_model.negThtSlowSteps[cobra_state.cobraidx[mask], :]

        assert ~np.any(np.isnan(theta_steps))

        return theta_steps, theta_offset

    def __get_phi_model(self, cobra_state, phi_delta, phi_fast):
        """
        Return the calibration model of the phi motor for each (virtual) cobra
        depending on the direction and speed.
        """

        # This is used at at least two locations in the original code:
        # - CobraCoach.calculateSteps
        # - calculus.interpolateMoves

        calib_steps = self.__calib_model.phiOffsets.shape[1]
        batch_shape = np.shape(cobra_state.cobraidx)
        shape = batch_shape + (calib_steps,)

        phi_offset = self.__calib_model.phiOffsets[cobra_state.cobraidx]
        phi_steps = np.full(shape, np.nan)

        # Depending on the direction of the movement and speed, select the correct model

        mask = (phi_delta >= 0) & phi_fast
        phi_steps[mask] = self.__calib_model.posPhiSteps[cobra_state.cobraidx[mask], :]

        mask = (phi_delta >= 0) & ~phi_fast
        phi_steps[mask] = self.__calib_model.posPhiSlowSteps[cobra_state.cobraidx[mask], :]

        mask = (phi_delta < 0) & phi_fast
        phi_steps[mask] = self.__calib_model.negPhiSteps[cobra_state.cobraidx[mask], :]

        mask = (phi_delta < 0) & ~phi_fast
        phi_steps[mask] = self.__calib_model.negPhiSlowSteps[cobra_state.cobraidx[mask], :]

        assert ~np.any(np.isnan(phi_steps))
        
        return phi_steps, phi_offset
    
    def __get_step_delay(self, cobra_state, theta_steps, phi_steps, max_steps):
        """
        Calculate the delay for the theta and phi motors in terms of number of steps.

        Parameters
        ----------
        cobra_state : SimpleNamespace
            Current cobra state.
        theta_steps : np.ndarray of int
            Number of steps for theta motor.
        phi_steps : np.ndarray of int
            Number of steps for phi motor.
        max_steps : int
            Maximum number of steps allowed.

        Returns
        -------
        theta_delay : np.ndarray of int
            Delay for the theta motor.
        phi_delay : np.ndarray of int
            Delay for the phi motor
        """

        # Set delay parameters for safer operation, same logic as in pfi.py

        theta_delay = np.zeros_like(theta_steps)
        phi_delay = np.zeros_like(phi_steps)
        
        mask_1 = (phi_steps > 0) & (theta_steps < 0)
        theta_delay[mask_1] = max_steps + theta_steps[mask_1]
        phi_delay[mask_1] = max_steps - phi_steps[mask_1]

        mask_2 = ~mask_1 & (phi_steps > 0) & (theta_steps > 0)
        theta_delay[mask_2] = 0
        phi_delay[mask_2] = max_steps - phi_steps[mask_2]

        mask_3 = ~mask_1 & ~mask_2 & (phi_steps < 0) & (theta_steps < 0)
        theta_delay[mask_3] = max_steps + theta_steps[mask_3]
        phi_delay[mask_3] = 0

        mask_4 = ~mask_1 & ~mask_2 & ~mask_3
        theta_delay[mask_4] = 0
        phi_delay[mask_4] = 0

        return theta_delay, phi_delay
        
    def __calculate_steps(self, cobra_state, theta_delta, phi_delta, theta_from, phi_from, theta_fast, phi_fast, max_steps):
        """
        This function is very closely based on CobraCoach.calculateSteps but with vectorization support.
        """

        # Select the correct theta model for each (virtual) cobra
        theta_model, theta_offset = self.__get_theta_model(cobra_state, theta_delta, theta_fast)
        phi_model, phi_offset = self.__get_phi_model(cobra_state, phi_delta, phi_fast)

        # Calculate the number of steps
        # --------------------------------------------------------------

        # Calculate the total number of motor steps for the theta movement
        # TODO: replace this with batch-mode interpolation

        steps_range = self.__batch_interp(np.stack([theta_from, (theta_from + theta_delta)], axis=-1),
                                          theta_offset,
                                          theta_model)
                    
        steps_range[~np.isfinite(steps_range)] = max_steps
        
        assert np.all(np.isfinite(steps_range))

        theta_steps = np.rint(steps_range[..., 1] - steps_range[..., 0]).astype(int)
        theta_from_steps = steps_range[..., 0]

        # --------------------------------------------------------------

        # Calculate the total number of motor steps for the phi movement
        steps_range = self.__batch_interp(np.stack([phi_from, (phi_from + phi_delta)], axis=-1),
                                          phi_offset,
                                          phi_model)
            
        steps_range[~np.isfinite(steps_range)] = max_steps

        assert np.all(np.isfinite(steps_range))

        phi_steps = np.rint(steps_range[..., 1] - steps_range[..., 0]).astype(int)
        phi_from_steps = steps_range[..., 0]

        # Update phi motor steps away from center in the beginning
        phi_steps_mask = phi_steps > 0
        phi_safe = cobra_state.phi_threshold - self.__calib_model.phiIn[cobra_state.cobraidx][phi_steps_mask] - np.pi
        phi_safe_mask = phi_from[phi_steps_mask] < phi_safe

        steps_range[phi_steps_mask][phi_safe_mask] = \
                self.__batch_interp(np.stack([phi_from[phi_steps_mask][phi_safe_mask], phi_safe[phi_safe_mask]], axis=-1),
                                    self.__calib_model.phiOffsets[cobra_state.cobraidx[phi_steps_mask][phi_safe_mask]],
                                    phi_model[phi_steps_mask][phi_safe_mask])

        assert np.all(np.isfinite(steps_range))

        safe_phi_steps = np.rint(steps_range[..., 1] - steps_range[..., 0]).astype(int)
        safe_phi_steps = np.minimum(safe_phi_steps, phi_steps)
        safe_phi_steps = np.minimum(safe_phi_steps, max_steps)

        # --------------------------------------------------------------
        
        # Determine the total number of steps
        too_many_steps = (np.abs(theta_steps) > max_steps) | (abs(phi_steps) > max_steps)
        more_theta_steps = too_many_steps & (np.abs(theta_steps) >= np.abs(phi_steps))
        more_phi_steps = ~more_theta_steps

        #

        mask_3 = (phi_steps > 0) & (theta_steps > 0)
        phi_steps[more_theta_steps & mask_3] = np.maximum(phi_steps - theta_steps + max_steps, safe_phi_steps)[more_theta_steps & mask_3]
        theta_steps[more_theta_steps & mask_3] = max_steps

        theta_steps[more_phi_steps & mask_3] = np.minimum(theta_steps, max_steps)[more_phi_steps & mask_3]
        phi_steps[more_phi_steps & mask_3] = max_steps

        #

        mask_4 = (phi_steps > 0) & (theta_steps < 0)
        phi_steps[more_theta_steps & mask_4] = np.maximum(phi_steps + theta_steps + max_steps, safe_phi_steps)[more_theta_steps & mask_4]
        theta_steps[more_theta_steps & mask_4] = -max_steps

        theta_steps[more_phi_steps & mask_4] = np.minimum(theta_steps + phi_steps - max_steps, 0)[more_phi_steps & mask_4]
        phi_steps[more_phi_steps & mask_4] = max_steps

        #

        mask_5 = (phi_steps < 0) & (theta_steps > 0)
        phi_steps[more_theta_steps & mask_5] = np.maximum(phi_steps, -max_steps)[more_theta_steps & mask_5]
        theta_steps[more_theta_steps & mask_5] = max_steps

        theta_steps[more_phi_steps & mask_5] = np.minimum(theta_steps, max_steps)[more_phi_steps & mask_5]
        phi_steps[more_phi_steps & mask_5] = -max_steps

        #

        mask_6 = (phi_steps < 0) & (theta_steps < 0)
        phi_steps[more_theta_steps & mask_6] = np.maximum(phi_steps, -max_steps)[more_theta_steps & mask_6]
        theta_steps[more_theta_steps & mask_6] = -max_steps

        theta_steps[more_phi_steps & mask_6] = np.minimum(theta_steps - phi_steps - max_steps, 0)[more_phi_steps & mask_6]
        phi_steps[more_phi_steps & mask_6] = -max_steps

        # The last two cases are different for phi and theta!

        mask_7 = (theta_steps > 0)
        theta_steps[more_theta_steps & mask_7] = max_steps

        # TODO: this might be a bit simpler, actually
        mask_8 = ~mask_3 & ~mask_4 & ~mask_5 & ~mask_6 & ~mask_7
        theta_steps[more_theta_steps & mask_8] = -max_steps

        #

        mask_7 = (phi_steps > 0)
        phi_steps[more_phi_steps & mask_7] = max_steps

        # TODO: this might be a bit simpler, actually
        mask_8 = ~mask_3 & ~mask_4 & ~mask_5 & ~mask_6 & ~mask_7
        phi_steps[more_phi_steps & mask_8] = -max_steps

        # --------------------------------------------------------------

        theta_to = self.__batch_interp(theta_from_steps + theta_steps, theta_model, theta_offset)
        phi_to = self.__batch_interp(phi_from_steps + phi_steps, phi_model, phi_offset)

        return theta_steps, phi_steps, theta_to - theta_from, phi_to - phi_from

    def __batch_interp(self, x, xp, fp):
        """
        Interpolates values in batch mode.
        """

        batch_shape = xp.shape[:-1]
        value_shape = x.shape[len(batch_shape):]

        r = np.empty(batch_shape + value_shape, dtype=float)
        for ix in np.ndindex(batch_shape):
            r[ix] = np.interp(x[ix], xp[ix], fp[ix])

        return r

    def __calculate_segments(self, angle_to, angle_from, mm, max_steps):
        """
        This is a vectorized version of cobra coach's calNSegments function.
        """

        moved = np.zeros_like(angle_to)
        n = np.zeros_like(angle_to, dtype=int)

        raise NotImplementedError()

        # if angle > 0:
        #     while moved < angle_to:
        #         idx = np.nanargmin(abs(mm[0]['angle'] - angle_from - moved))
        #         moved += mm[0,idx]['speed'] * maxSteps
        #         n += 1
        # elif angle < 0:
        #     while moved > angle:
        #         # Asking program to look for negtive speed in on-time map. Otherwise, it will fall into
        #         # infinite loop.
        #         #idx=np.nanargmin(abs(mm[1,np.where(mm[1]['speed']<0)]['angle']- fromAngle - moved))
        #         idx = np.nanargmin(abs(mm[1]['angle'] - fromAngle - moved))
        #         moved += mm[1,idx]['speed'] * maxSteps
        #         n += 1

        # return n

    def find_associations(self, *coords, mask=None):
        
        # TODO: this is deprecated, not using netflow
        #       move elsewhere or delete

        ctype, coords = normalize_coords(*coords)
        mask = mask if mask is not None else ()

        xy, fov_mask = self.__projection.world_to_pixel(coords[mask])
    
        b = self.__bench
        centers = np.array([b.cobras.centers.real, b.cobras.centers.imag]).T

        if fov_mask.sum() > 10:
            # Use kD-tree to find targets within R_max and
            # filter out results that are within R_min
            tree = KDTree(xy[fov_mask], leaf_size=2)
            outer = tree.query_radius(centers, b.cobras.rMax)
            inner = tree.query_radius(centers, b.cobras.rMin)

            assoc = [ np.setdiff1d(i, o, assume_unique=True) for i, o in zip(outer, inner) ]

            # Assoc is now a list with size equal to the number of cobras.
            # Each item in an integer array with the list of matching indexes into
            # the coords[mask][fov_mask] array.

            # Update associations to index into coords[mask] instead of coords[mask][fov_mask]
            index = np.arange(coords[mask].shape[0], dtype=np.int64)[fov_mask]
            assoc = [ index[a] for a in assoc ]
        else:
            # Do a 1-1 matching instead?
            assoc = [ np.array([], dtype=np.int64) for i in range(len(b.cobras.centers))]
        
        return Associations(assoc)
    
    def __plot_corners(self, ax):
        pass

    def __get_outline(self, ids, res, native_frame=None, projection=None):
        
        projection = projection if projection is not None else self.__projection

        # Calculate corners
        xy = []
        for i in ids:
            xy.append((self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag))
        xy.append(xy[0])    # wrap up
        xy = np.array(xy)

        # Calculate arcs
        axy = []
        t = np.linspace(0, 1, res)
        for xy0, xy1 in zip(xy[:-1], xy[1:]):
            axy.append(xy0 + t[:, np.newaxis] * (xy1 - xy0)[np.newaxis, :])
        
        axy = np.concatenate(axy)
        if native_frame == 'world':
            coords, mask = projection.pixel_to_world(axy)
        elif native_frame == 'pixel':
            coords = axy
            mask = np.full(axy.shape[:-1], True, dtype=bool)
        
        return coords, mask

    def plot_focal_plane(self, ax: plt.Axes, diagram, res=None, projection=None, **kwargs):
        corners = kwargs.pop('corners', False)
        blocks = kwargs.pop('blocks', False)
        fill = kwargs.pop('fill', False)

        res = res if res is not None else 36
        native_frame = diagram._get_native_frame()

        if not fill:
            style = styles.solid_line(**kwargs)
            def plot_outline(ids):
                for ii in ids:
                    xy, mask = self.__get_outline(ii, res, native_frame=native_frame, projection=projection)
                    diagram.plot(ax, xy, mask=mask, native_frame=native_frame,**style)

            if corners:
                plot_outline(self.CORNERS)
            
            if blocks:
                plot_outline(self.BLOCKS)
        else:
            style = styles.red_fill(**kwargs)
            for ii in self.CORNERS:
                xy, mask = self.__get_outline(ii, res, native_frame=native_frame, projection=projection)
                diagram.fill(ax, xy, mask=mask, native_frame=native_frame, **style)

    def plot_cobras(self, ax, diagram, data=None, cmap='viridis', vmin=None, vmax=None,
                    scalex=True, scaley=True,
                    **kwargs):
        """
        Plots the fibers color-coded by the data vector
        """

        cmap = get_cmap(cmap) if not isinstance(cmap, Colormap) else cmap
        native_frame = diagram._get_native_frame()

        if data is not None:
            style = styles.closed_circle(**kwargs)
            vmin, vmax = find_plot_limits(data, vmin=vmin, vmax=vmax)
            color = get_plot_normalized_color(cmap, data, vmin=vmin, vmax=vmax, alpha=0.5)
        else:
            style = styles.open_circle(**kwargs)
            vmin, vmax = None, None
            color = None

        # Draw the cobra patrol areas with actual circles - this is slow
        if False:
            for i in range(len(self.__bench.cobras.centers)):
                # TODO: we are cheating here and assume that the patrol region will appear as a cirle
                #       in the plots, regardless of the actual projection

                # TODO: This projects cobras into FP coordinates so it won't work 
                #       when we draw a FOV plot with some projection defines on the axis

                # TODO: review plotting broken cobras

                if isinstance(color, Iterable):
                    style['facecolor'] = color[i]
                elif color is not None:
                    style['facecolor'] = color

                xy = (self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag)
                rr = (self.__bench.cobras.centers[i].real + self.__bench.cobras.rMax[i], self.__bench.cobras.centers[i].imag)
                xy = diagram.project_coords(ax, xy, native_frame='pixel')
                rr = diagram.project_coords(ax, rr, native_frame='pixel')
                r = rr[0] - xy[0]

                c = Circle(xy, r, **style)
                diagram.add_patch(ax, c)

        # Draw the cobras with symbols
        if True:
            ax.scatter(self.__bench.cobras.centers.real,
                        self.__bench.cobras.centers.imag,
                        marker='o',
                        s=8,
                        c=color, **style)

        return ScalarMappable(cmap=cmap, norm=Normalize(vmin, vmax))

    # TODO: plot black dots

    def plot_fiber_numbers(self, ax, diagram, **kwargs):

        # TODO: figure out what coordinates and projections to use

        for i in range(len(self.__bench.cobras.centers)):
            x, y = (self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag)
            x, y = diagram.project_coords(ax, x, y, native_frame='pixel')
            ax.text(x, y, str(i), ha='center', va='center', **kwargs)

    def plot_corners(self, ax, **kwargs):
        pass

    def generate_cobra_location_labels(self, ntheta=6):
        return self.__generate_cobra_location_labels_voronoi(ntheta=ntheta)

    def __generate_cobra_location_labels_voronoi(self, ntheta=None):

        # Get the focal plane coordinates of the cobras
        uv = np.stack([self.__bench.cobras.centers.real, self.__bench.cobras.centers.imag], axis=-1)

        # Put down some points in the focal plane and flag the cobras with the index of the closest point

        points = []
        points.append([[0, 0]])

        phi = np.radians(np.linspace(0, 360, 6, endpoint=False) + 15 )
        points.append(100 * np.stack([np.cos(phi), np.sin(phi)], axis=-1))

        phi = np.radians(np.linspace(0, 360, 13, endpoint=False) - 15)
        points.append(175 * np.stack([np.cos(phi), np.sin(phi)], axis=-1))

        xy = np.concatenate(points, axis=0)

        d = distance_matrix(xy, uv)

        labels = np.argmin(d, axis=0)

        return labels
    
    def __generate_cobra_location_labels_sections(self, ntheta=6):
        """
        Generate integer labels for each cobra that organize the cobras into groups that
        are uniform in the focal plane.
        
        Parameters
        ----------
        ntheta : int
            The number of groups to divide the circle into.
        """

        # Get the focal plane coordinates of the cobras
        x, y = self.__bench.cobras.centers[:].real, self.__bench.cobras.centers[:].imag

        # Convert to polar coordinates around the center of the focal plane
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)

        # Assign labels to the cobras based on the polar coordinate
        theta_bins = np.linspace(-np.pi, np.pi, ntheta + 1)
        r_bins = np.array([0, 150, 240])

        theta_labels = np.digitize(theta, theta_bins, right=False) - 1
        r_labels = np.digitize(r, r_bins, right=False) - 1
        cobra_location_labels = (r_bins.size - 1) * theta_labels + r_labels

        # Add one more label in the center
        cobra_location_labels[r < 60] = cobra_location_labels.max() + 1

        return cobra_location_labels
    
    def generate_cobra_instrument_labels(self, ngroups=8):
        """
        Generate integer labels for each cobra that organize the cobras into groups that
        are uniform along the slit in each spectrograph.
        
        Parameters
        ----------
        ngroups : int
            The number of groups to divide the cobras into, for each spectrograph.
        """

        ncobras = len(self.__bench.cobras.centers)
        mask = self.__fiber_map.cobraId != FiberIds.MISSING_VALUE       # Non-engineering fibers
        
        cobra_instrument_labels = np.zeros(ncobras, dtype=int)
        for s in np.arange(np.unique(self.__fiber_map.spectrographId).max()):
            mask_s = mask & (self.__fiber_map.spectrographId == s + 1)     # Fibers connected to spectrograph s + 1
            fiber_count = np.sum(mask_s)                                # Number of fibers connected to spectrograph s + 1           
            cobra_instrument_labels[self.__fiber_map.cobraId[mask_s] - 1] = s * ngroups \
                + np.floor(np.arange(fiber_count) / fiber_count * ngroups).astype(int)

        return cobra_instrument_labels    

