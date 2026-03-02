import numpy as np
import astropy.units as u
from datetime import datetime, timedelta, tzinfo

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ... import Isochrone
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection, IsochroneSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ...config.pmap import PMapConfig
from ...config.sample import SampleConfig
from ..ids import *
from .gc import GC

from ...setup_logger import logger

class NGC1904(GC):
    def __init__(self):
        ID = 'NGC1904'
        name = 'NGC1904 (M79)'

        pos = [ 210.025, 14.5 ] * u.deg                       # Evan
        rad = 120 * u.arcmin
        DM, DM_err = 19.11, 0.008                             # Oakes et al. (2022)
        pm = [ -0.387, -1.064 ] * u.mas / u.yr                # Evan
        pm_err = [ 0.122, 0.098 ] * u.mas / u.yr
        RV, RV_err = (99.0, 2.1) * u.kilometer / u.second     # Simbad

        # # SSP pointings
        # ra0 = [ 210.5, 209.6, 210.1, 210.1 ] * u.deg
        # dec0 = [ 14.5, 14.5, 14.15, 14.8 ] * u.deg
        # pa0 = [ 30, 30, 30, 30 ] * u.deg

        # # Engineering run
        # # ra0 = [ 210.025 ] * u.deg
        # # dec0 = [ 14.5 ] * u.deg
        # # pa0 = [ 30 ] * u.deg

        # pointings = {
        #     SubaruPFI: [ Pointing((ra, dec), posang=pa, stage=0)
        #                  for ra, dec, pa in zip(ra0, dec0, pa0) ]
        # }

        pointings = {
            SubaruPFI: [
                Pointing.from_relative_pos(pos, sep=0.22, dir=40, posang=0, stage=0, priority=1),
                Pointing.from_relative_pos(pos, sep=-0.22, dir=40, posang=0, stage=0, priority=1),
                Pointing(pos, posang=0, stage=1, priority=4),
            ]
        }

        super().__init__(ID, name, ID_PREFIX_BOOTES,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)

