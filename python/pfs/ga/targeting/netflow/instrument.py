import os

import pfs.utils
from pfs.utils.coordinates import DistortionCoefficients as DCoeff
from pfs.utils.coordinates.CoordTransp import CoordinateTransform
import pfs.instdata
from ics.cobraOps.Bench import Bench
from pfs.utils.fiberids import FiberIds
from ics.cobraCharmer.pfiDesign import PFIDesign
from ics.cobraOps.BlackDotsCalibrationProduct import BlackDotsCalibrationProduct
from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach

from .setup_logger import logger
from .util import *

class Instrument():
    """
    Provides function to initialize the PFS instrument configurations.

    Configuration steps are adopted from the PFS design tool:
    ets_pointing/src/pfs_design_tool/pointing_utils/nfutils.py

    instrument_options: dict
        A dictionary with the following keys:
        - fiberids_path: str
            Path to the grand fiber map file.
        - instdata_path: str
            Path to the instrument configuration.
        - blackdots_path: str
            Path to the black dots calibration file.
        
    """

    # Path to grand fiber map file
    DEFAULT_FIBERIDS_PATH = os.path.join(os.path.dirname(pfs.utils.__file__), '../../../data/fiberids')
    
    # Path to instrument configuration
    DEFAULT_INSTDATA_PATH = os.path.join(os.path.dirname(pfs.instdata.__file__), '../../../')

    # Path to black dots calibration file
    DEFAULT_BLACK_DOTS_PATH = os.path.join(DEFAULT_INSTDATA_PATH, "data/pfi/dot", "black_dots_mm.csv")

    # Path to the cobra coach directory
    DEFAULT_COBRA_COACH_DIR = os.path.expandvars('/tmp/$USER/cobraCoach')

    def __init__(self, instrument_options=None):
        self.__instrument_options = camel_to_snake(instrument_options)

    def __get_instrument_option(self, key, default=None):
        """
        Return an option from the instrument_options dictionary, if it exists,
        otherwise return the default value.
        """

        if self.__instrument_options is not None and key in self.__instrument_options:
            return self.__instrument_options[key]
        else:
            return default

    def __get_grand_fiber_map(self):
        fiberids_path = self.__get_instrument_option('fiberids_path', Instrument.DEFAULT_FIBERIDS_PATH)
        fibermap = FiberIds(path=fiberids_path)
        return fibermap

    def __get_default_bench(self):
        return Bench(layout='full')

    def __get_configured_bench(self):
        """
        Create a bench object with real instrument configuration.
        """
        
        pfs_instdata_path = self.__get_instrument_option('instdata_path', Instrument.DEFAULT_INSTDATA_PATH)
        pfs_black_dots_path = self.__get_instrument_option('blackdots_path', Instrument.DEFAULT_BLACK_DOTS_PATH)
        cobra_coach_dir = self.__get_instrument_option('cobra_coach_dir', Instrument.DEFAULT_COBRA_COACH_DIR)
        cobra_coach_module_version = self.__get_instrument_option('cobra_coach_module_version', None)
        spectgrograph_modules = self.__get_instrument_option('spectrograph_modules', [1, 2, 3, 4])
        black_dot_radius_margin = self.__get_instrument_option('black_dot_radius_margin', 1.0)

        # Create the cobra coach temp directory if it does not exist
        if not os.path.isdir(cobra_coach_dir):
            os.makedirs(cobra_coach_dir, exist_ok=True)
            logger.info(f"Created cobra coach temp directory: {cobra_coach_dir}")

        # This might be required somewhere in the libraries
        os.environ["PFS_INSTDATA_DIR"] = pfs_instdata_path
        
        cobraCoach = CobraCoach("fpga", loadModel=False, trajectoryMode=True, rootDir=cobra_coach_dir)
        cobraCoach.loadModel(version="ALL", moduleVersion=cobra_coach_module_version)
        calibrationProduct = cobraCoach.calibModel

        # Set some dummy center positions and phi angles for those cobras that have
        # zero centers
        zeroCenters = calibrationProduct.centers == 0
        calibrationProduct.centers[zeroCenters] = np.arange(np.sum(zeroCenters)) * 300j
        calibrationProduct.phiIn[zeroCenters] = -np.pi
        calibrationProduct.phiOut[zeroCenters] = 0
        logger.info("Cobras with zero centers: %i" % np.sum(zeroCenters))

        # Use the median value link lengths in those cobras with zero link lengths
        zeroLinkLengths = (calibrationProduct.L1 == 0) | (calibrationProduct.L2 == 0)
        calibrationProduct.L1[zeroLinkLengths] = np.median(calibrationProduct.L1[~zeroLinkLengths])
        calibrationProduct.L2[zeroLinkLengths] = np.median(calibrationProduct.L2[~zeroLinkLengths])
        logger.info("Cobras with zero link lengths: %i" % np.sum(zeroLinkLengths))

        # Use the median value link lengths in those cobras with too long link lengths
        tooLongLinkLengths = (calibrationProduct.L1 > 100) | (calibrationProduct.L2 > 100)
        calibrationProduct.L1[tooLongLinkLengths] = np.median(calibrationProduct.L1[~tooLongLinkLengths])
        calibrationProduct.L2[tooLongLinkLengths] = np.median(calibrationProduct.L2[~tooLongLinkLengths])
        logger.info("Cobras with too long link lengths: %i" % np.sum(tooLongLinkLengths))

        # Limit spectral modules
        gfm = self.__get_grand_fiber_map()
        cobra_ids_use = np.array([], dtype=np.uint16)
        for sm in spectgrograph_modules:
            cobra_ids_use = np.append(cobra_ids_use, gfm.cobrasForSpectrograph(sm))

        # Set Bad Cobra status for unused spectral modules
        for cobra_id in range(calibrationProduct.nCobras):
            if cobra_id not in cobra_ids_use:
                calibrationProduct.status[cobra_id] = ~PFIDesign.COBRA_OK_MASK
        
        # Get the black dots calibration product
        blackDotsCalibrationProduct = BlackDotsCalibrationProduct(pfs_black_dots_path)

        bench = Bench(
            layout="calibration",                       # Use the layout from the calibration product
            calibrationProduct=calibrationProduct,
            blackDotsCalibrationProduct=blackDotsCalibrationProduct,
            blackDotsMargin=black_dot_radius_margin,
        )

        return bench
    
    def get_fiber_map(self):
        return self.__get_grand_fiber_map()
    
    def get_bench(self):
        layout = self.__get_instrument_option('layout', 'full')

        if layout == 'full':
            return self.__get_default_bench()
        elif layout == 'calibration':
            return self.__get_configured_bench()

    def radec_to_altaz(self,
                       ra, dec, posang, obs_time,
                       epoch=2000.0, pmra=0.0, pmdec=0.0, parallax=1e-6, rv=0.0):

        az, el, inr = DCoeff.radec_to_subaru(ra, dec, posang, obs_time,
                                             epoch=epoch, pmra=pmra, pmdec=pmdec, par=1e-6)

        return az, el, inr
    
    def radec_to_fp_pos(self,
                        pointing,
                        ra, dec,
                        epoch=2000.0, pmra=None, pmdec=None, parallax=None, rv=None):
        
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