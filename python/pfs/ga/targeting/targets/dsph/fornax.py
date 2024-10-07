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

class Fornax(DSphGalaxy):
    def __init__(self):
        ID = 'for'
        pos = [ 39.9971, -34.4492 ] * u.deg                     # Evan
        rad = np.nan * u.arcmin
        DM, DM_err = 20.77, 0.05                                # Oakes et al. (2022)
        pm = [ 0.381, -0.358 ] * u.mas / u.yr                   # Evan
        pm_err = [ 0.001, 0.002 ] * u.mas / u.yr
        RV, RV_err = (-291.0, 0.1) * u.kilometer / u.second     # Simbad

        ra0 = [ 39.5, 40.4, 40.3, 39.5, 40.9, 40.5, 39.4, 39.1 ] * u.deg
        dec0 = [ -34.2, -34.1, -34.8, -34.9, -33.8, -35.1, -34.0, -35.2 ] * u.deg
        pa0 = [ 30, 30, 30, 30, 30, 30, 30, 30 ] * u.deg

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