import astropy.units as u
from datetime import datetime, timedelta, tzinfo
from collections import defaultdict

from sklearn import logger

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ...config.pmap import PMapConfig
from ...config.sample import SampleConfig
from ..galaxy import Galaxy
from ..ids import *

from ...setup_logger import logger

from .m31 import M31

class M33(M31):

    FIELDS = [
        # group_name, field_name,  ra,             dec, posang, stage, priority

        # Fields prepared for 2026-11

        ('sector_0', 'm33_center',   '01:33:49.99',  '+30:39:37.1', 50.0, 0, 0),
        ('sector_0', 'm33_SW',       '01:28:49.72',  '+30:15:39.5', 50.0, 0, 0),
        ('sector_0', 'm33_NW',       '01:29:42.30',  '+31:23:43.5', 50.0, 0, 0),
        ('sector_0', 'm33_N',        '01:34:46.38',  '+31:47:33.4', 50.0, 0, 0),
        ('sector_0', 'm33_NE',       '01:38:52.72',  '+31:02:51.1', 50.0, 0, 0),
        ('sector_0', 'm33_SE',       '01:37:53.92',  '+29:55:01.7', 50.0, 0, 0),
        ('sector_0', 'm33_S',        '01:32:54.91',  '+29:31:39.2', 50.0, 0, 0),
        
    ]

    def __init__(self, sector=None, field=None):
    
        ID = 'm33'
        name = 'M33'
        self.sector = sector
        self.field = field

        pos = [ '01h 33m 50.90s', '+30d 39m 36.6s' ]
        rad = 3 * u.deg
        DM, DM_err = 24.67, 0.06   #Savino et al. (2022, ApJ, 938, 101)
        pm = [ 0.0511, 0.0114 ] * u.mas / u.yr    #Rusterucci et al. (2024, A&A, 692, A30)
        pm_err = [ 0.0568, 0.0500 ] * u.mas / u.yr
        RV, RV_err = (-179, 1) * u.kilometer / u.second  #Koch et al. (2018, MNRAS, 479, 2505)

        pointings_by_sector = defaultdict(dict)
        for sc, fl, ra, dec, pa0, st, pri in self.FIELDS:
            pointings_by_sector[sc][fl] = Pointing(
                Angle(ra, unit=u.hourangle),
                Angle(dec, unit=u.deg),
                posang=pa0, priority=pri, stage=st,
                exp_time=30*60, nvisits=10,
                label=fl)
        
        if sector is None and field is None:
            # All pointings
            pointings = {
                SubaruPFI: [ pointing
                             for s in pointings_by_sector
                             for f, pointing in pointings_by_sector[s].items() ]
            }
        elif sector is not None and field is None:
            # Use a single sector with all its pointings
            pointings = {
                SubaruPFI: [ pointing
                             for f, pointing in pointings_by_sector[sector].items() ]
            }
        else:
            # Use a single pointing
            pointings = {
                SubaruPFI: [ pointings_by_sector[sector][field] ]
            }

        # Skip M31.__init__ (different signature and M31 constants) and
        # initialize the Galaxy base class directly
        Galaxy.__init__(self, ID, name, ID_PREFIX_M33,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        # CMD and CCD definitions with color and magnitude limits for M33

        hsc = SubaruHSC.photometry()
        self._hsc_cmd = CMD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
            MagnitudeAxis(hsc.magnitudes['g'], limits=(15.5, 24.5))
        ])
        self._hsc_ccd = CCD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
            ColorAxis( Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), limits=(-0.5, 0.5))
        ])

        gaia = Gaia.photometry()
        self._gaia_cmd = CMD([
            ColorAxis(Color([gaia.magnitudes['bp'], gaia.magnitudes['rp']]), limits=(0, 3)),
            MagnitudeAxis(gaia.magnitudes['g'], limits=(11, 22))
        ])