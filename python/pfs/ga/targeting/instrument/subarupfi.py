import os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.cm import get_cmap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from sklearn.neighbors import KDTree

import pfs.utils
from pfs.utils.coordinates import DistortionCoefficients as DCoeff
from pfs.utils.coordinates.CoordTransp import CoordinateTransform
import pfs.instdata
from ics.cobraOps.Bench import Bench
from pfs.utils.fiberids import FiberIds
from pfs.utils.butler import Butler
from ..external import BlackDotsCalibrationProduct
from ..external import PFIDesign
from ..external import CobraCoach

from ..setup_logger import logger
from ..util import *
from ..util.config import *
from ..diagram import styles
from ..allocation.associations import Associations
from ..allocation.fiberallocator import FiberAllocator
from ..projection import Pointing
from .instrument import Instrument
from .subaruwfc import SubaruWFC
from .cobraangleflags import CobraAngleFlags

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

    def __init__(self, projection=None, instrument_options=None, orig=None):
        from ..config.instrumentoptionsconfig import InstrumentOptionsConfig
        
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
        self.__bench, self.__cobra_coach, self.__calib_model = self.__create_bench()

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
        fiber_map = FiberIds(path=fiberids_path)

        logger.info(f"Loading the grand fiber map from `{os.path.abspath(fiberids_path)}`")

        # Load the list of blocked fibers
        config_root = os.path.join(self.__get_instrument_option(self.__instrument_options.instdata_path, SubaruPFI.DEFAULT_INSTDATA_PATH), 'data')
        butler = Butler(configRoot=config_root)
        blocked_fibers = butler.get('fiberBlocked').set_index('fiberId')

        logger.info(f"Loaded configuration for {len(fiber_map.fiberId)} fibers.")
        logger.info(f"Number of fibers connected to cobras: {(fiber_map.cobraId != fiber_map.MISSING_VALUE).sum()}")
        logger.info(f"Number of science fibers: {(fiber_map.scienceFiberId < fiber_map.ENGINEERING).sum()}")

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
    
    def __create_bench(self):
        layout = self.__get_instrument_option(self.__instrument_options.layout, 'full')

        if layout == 'full':
            return self.__create_default_bench()
        elif layout == 'calibration':
            return self.__create_configured_bench()
        else:
            raise NotImplementedError()

    def __create_default_bench(self):
        """
        Create a bench object with default instrument configuration. This is
        not the real configuration and the PFI is rotated by 30 degrees.
        """
        return Bench(layout='full'), None, None
    
    def __create_configured_bench(self):
        """
        Create a bench object with real instrument configuration.
        """
        
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
        
        # Get the black dots calibration product
        black_dots_calibration_product = BlackDotsCalibrationProduct(pfs_black_dots_path)

        bench = Bench(
            layout="calibration",                       # Use the layout from the calibration product
            calibrationProduct=calib_model,
            blackDotsCalibrationProduct=black_dots_calibration_product,
            blackDotsMargin=black_dot_radius_margin,
        )

        return bench, cobra_coach, calib_model
    
    def __fix_calib_model(self, calib_model):
        """
        Implemnents some calib model fixes that are not required anymore.
        """

        ignore_calibration_errors = self.__get_instrument_option(self.__instrument_options.ignore_calibration_errors, False)

        # Set some dummy center positions and phi angles for those cobras that have
        # zero centers
        zero_centers = calib_model.centers == 0
        if zero_centers.any():
            msg = f"Cobras with zero centers: {np.sum(zero_centers)}"
            if not ignore_calibration_errors:
                raise ValueError(msg)
            else:
                calib_model.centers[zero_centers] = np.arange(np.sum(zero_centers)) * 300j
                calib_model.phiIn[zero_centers] = -np.pi
                calib_model.phiOut[zero_centers] = 0
                logger.warning(msg)
                
        # Use the median value link lengths in those cobras with zero link lengths
        zero_link_lengths = (calib_model.L1 == 0) | (calib_model.L2 == 0)
        if zero_link_lengths.any():
            msg = f"Cobras with zero link lengths: {np.sum(zero_link_lengths)}"
            if not ignore_calibration_errors:
                raise ValueError(msg)
            else:
                calib_model.L1[zero_link_lengths] = np.median(calib_model.L1[~zero_link_lengths])
                calib_model.L2[zero_link_lengths] = np.median(calib_model.L2[~zero_link_lengths])
                logger.warning(msg)

        # Use the median value link lengths in those cobras with too long link lengths
        too_long_link_lengths = (calib_model.L1 > 100) | (calib_model.L2 > 100)
        if too_long_link_lengths.any():
            msg = f"Cobras with too long link lengths: {np.sum(too_long_link_lengths)}"
            if not ignore_calibration_errors:
                raise ValueError(msg)
            else:
                calib_model.L1[too_long_link_lengths] = np.median(calib_model.L1[~too_long_link_lengths])
                calib_model.L2[too_long_link_lengths] = np.median(calib_model.L2[~too_long_link_lengths])
                logger.warning(msg)

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
    
    def radec_to_fp_pos(self,
                        ra, dec,
                        pointing,
                        epoch=2000.0, pmra=None, pmdec=None, parallax=None, rv=None):
        
        """
        Convert celestial coordinate to focal plane position.
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
    
    def __get_reshaped_cobra_config(self, ndim, cobraidx):
        """
        Reshape the cobra configuration to match the shape of the input data.
        """
        
        batch_shape = (Ellipsis,) + (ndim - cobraidx.ndim) * (None,)
        centers = self.__bench.cobras.centers[cobraidx][batch_shape]
        bad_cobra = self.__bench.cobras.hasProblem[cobraidx][batch_shape]
        L1 = self.__bench.cobras.L1[cobraidx][batch_shape]
        L2 = self.__bench.cobras.L2[cobraidx][batch_shape]

        return batch_shape, bad_cobra, centers, L1, L2
    
    def fp_pos_to_rel_fp_pos(self, fp_pos, cobraidx):
        """
        Convert focal plane positions to relative focal plane positions. It also returns
        a few precomputed parameters.
        """

        batch_shape, bad_cobra, centers, L1, L2 = self.__get_reshaped_cobra_config(fp_pos.ndim, cobraidx)

        rel_fp_pos = fp_pos - centers
        d = np.abs(rel_fp_pos)
        d_2 = d * d

        return batch_shape, bad_cobra, rel_fp_pos, centers, L1, L2, d, d_2 
        
    def fp_pos_to_cobra_angles(self, fp_pos, cobraidx):
        """
        Convert focal plane positions to cobra angles.

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
        flags : np.ndarray
            Flags with shape (N, 2) where N is the number of targets indicating the validity of the solutions.
        """

        batch_shape, bad_cobra, rel_fp_pos, centers, L1, L2, d, d_2 = \
            self.fp_pos_to_rel_fp_pos(fp_pos, cobraidx)

        L1_2 = L1 ** 2
        L2_2 = L2 ** 2

        arp = np.angle(rel_fp_pos)

        # This might generate warnings because we don't filter on unreachable targets yet
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
        Convert cobra angles to focal plane position.
        """

        batch_shape, bad_cobra, centers, L1, L2 = self.__get_reshaped_cobra_config(theta.ndim, cobraidx)

        phi_in = self.__bench.cobras.phiIn[cobraidx][batch_shape]
        phi_out = self.__bench.cobras.phiOut[cobraidx][batch_shape]
        theta0 = self.__bench.cobras.tht0[cobraidx][batch_shape]
        theta1 = self.__bench.cobras.tht1[cobraidx[batch_shape]]

        # Do some sanity checks
        # TODO: move this to Netflow instead

        theta_range = (theta1 - theta0 + np.pi) % (2 * np.pi) + np.pi
        if np.any(0 > theta) or np.any(theta_range < theta):
            logger.error('Some theta angles are out of range')
        
        phi_range = phi_out - phi_in
        if np.any(0 > phi) or np.any(phi_range < phi):
            logger.error('Some phi angles are out of range')

        ang1 = theta + theta0
        ang2 = ang1 + phi + phi_in
        fp_pos = centers + L1 * np.exp(1j * ang1) + L2 * np.exp(1j * ang2)

        return fp_pos
    
    def simulate_trajectories(self, theta, phi, cobraidx, nsteps=100, two_step_mode=False):
        raise NotImplementedError()

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

        style = styles.solid_line(**kwargs)
        res = res if res is not None else 36
        native_frame = diagram._get_native_frame()

        def plot_outline(ids):
            for ii in ids:
                xy, mask = self.__get_outline(ii, res, native_frame=native_frame, projection=projection)
                diagram.plot(ax, xy, mask=mask, native_frame=native_frame, **style)

        if corners:
            plot_outline(self.CORNERS)
        
        if blocks:
            plot_outline(self.BLOCKS)           

    def plot_cobras(self, ax, diagram, data=None, cmap='viridis', vmin=None, vmax=None, **kwargs):
        """
        Plots the fibers color-coded by the data vector
        """

        cmap = get_cmap(cmap) if not isinstance(cmap, Colormap) else cmap
        native_frame = diagram._get_native_frame()

        if data is not None:
            style = styles.closed_circle(**kwargs)
            vmin, vmax = find_plot_limits(data, vmin=vmin, vmax=vmax)
            color = get_plot_normalized_color(cmap, data, vmin=vmin, vmax=vmax)
        else:
            style = styles.open_circle(**kwargs)
            vmin, vmax = None, None
            color = None

        for i in range(len(self.__bench.cobras.centers)):
            # TODO: we are cheating here and assume that the patrol region will appear as a cirle
            #       in the plots, regardless of the actual projection

            # TODO: This projects cobras into FP coordinates so it won't work 
            #       when we draw a FOV plot with some projection defines on the axis

            # TODO: review plotting broken cobras

            xy = (self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag)
            rr = (self.__bench.cobras.centers[i].real + self.__bench.cobras.rMax[i], self.__bench.cobras.centers[i].imag)
            xy = diagram.project_coords(ax, xy, native_frame='pixel')
            rr = diagram.project_coords(ax, rr, native_frame='pixel')
            r = rr[0] - xy[0]

            if isinstance(color, Iterable):
                style['facecolor'] = color[i]
            elif color is not None:
                style['facecolor'] = color

            c = Circle(xy, r, **style)
            ax.add_patch(c)

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
        """
        Generate integer labels for each cobra that organize the cobras into groups that
        are uniform in the focal plane.
        
        Parameters
        ----------
        ngroups : int
            The number of groups to divide the cobras into, for each spectrograph.
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

