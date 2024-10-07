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

class Draco(DSphGalaxy):
    def __init__(self):
        ID = 'dra'
        pos = [ 260.051666667, 57.9152777778 ] * u.deg
        rad = np.nan * u.arcmin
        DM, DM_err = 19.55, 0.0                                 # Wikipedia
        pm = [ 0.046, -0.188 ] * u.mas / u.yr                   # Evan
        pm_err = [ 0.006, 0.006 ] * u.mas / u.yr
        RV, RV_err = (-291.0, 0.1) * u.kilometer / u.second     # Simbad

        ra0 = [ 259.2, 260.9, 260.1, 260.0 ] * u.deg
        dec0 = [ 57.871, 57.957, 57.77, 58.05 ] * u.deg
        pa0 = [ 0, 0, 30, 30 ] * u.deg

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