import numpy as np
import astropy.units as u

from ...util.args import *
from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from ...photometry import Photometry, Magnitude, Color
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from .dsphgalaxy import DSphGalaxy

class NGC6822(DSphGalaxy):
    def __init__(self):
        ID = 'ngc6822'
        pos = [ 296.234, -14.7976 ] * u.deg
        rad = np.nan * u.arcmin
        DM, DM_err = np.nan, np.nan
        pm = [ np.nan, np.nan ] * u.mas / u.yr     
        pm_err = [ 0.122, 0.098 ] * u.mas / u.yr
        RV, RV_err = (np.nan, np.nan) * u.kilometer / u.second

        ra0 = [ 296.235 ] * u.deg
        dec0 = [ -14.789 ] * u.deg
        pa0 = [ 30 ] * u.deg

        pointings = {
            SubaruPFI: [ Pointing((ra, dec), posang=pa) for ra, dec, pa in zip(ra0, dec0, pa0) ]
        }

        super().__init__(ID,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        hsc = SubaruHSC.photometry()
        self.__hsc_cmd = CMD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
            MagnitudeAxis(hsc.magnitudes['g'], limits=(15.5, 24.5))
        ])
        self.__hsc_ccd = CCD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
            ColorAxis( Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), limits=(-0.5, 0.5))
        ])

        gaia = Gaia.photometry()
        self.__gaia_cmd = CMD([
            ColorAxis(Color([gaia.magnitudes['bp'], gaia.magnitudes['rp']]), limits=(0, 3)),
            MagnitudeAxis(gaia.magnitudes['g'], limits=(11, 22))
        ])