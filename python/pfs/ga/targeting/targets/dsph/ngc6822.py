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
        rad = 120. * u.arcmin
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

    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=18, faint=23.0):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self.__hsc_cmd
        ccd = self.__hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], None, 2.75).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= (
                #ColorSelection(ccd.axes[0], 0.12, 0.5).apply(catalog, observed=observed)

                #| ColorSelection(ccd.axes[1], 0.1, None).apply(catalog, observed=observed)
                #& ColorSelection(ccd.axes[0], None, 1.65).apply(catalog, observed=observed)
                
                #| LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed)
                ColorSelection(ccd.axes[0], None, 0.0).apply(catalog, observed=observed)

                | LinearSelection(ccd.axes, [0.25, 1.0], 0.25, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], 0.1, 1.5).apply(catalog, observed=observed)
                
                | LinearSelection(ccd.axes, [-0.33, 1.0], -0.62, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], 1.5, None).apply(catalog, observed=observed)
            )

        # Probability-based cut (map) - nonzero membership probability
        if probcut is not None:
            mask &= probcut.apply(catalog, observed=observed)

        # Allow blue
        if blue:
            mask |= (
                ColorSelection(ccd.axes[0], None, 0.12).apply(catalog, observed=observed)
            )

        # Always impose faint and bright magnitude cuts
        mask &= MagnitudeSelection(cmd.axes[1], bright, faint).apply(catalog, observed=observed)

        return mask

    def assign_priorities(self, catalog: Catalog, selection=None, mask=None):
        """Assign priority classes based on photometry"""

        mask = mask if mask is not None else np.s_[:]

        g0, i0, gi0 = self._get_hsc_dered_mags_colors(catalog, mask)
        clg = catalog.data['clg'][mask]
        cli = catalog.data['cli'][mask]

        priority = np.full(g0.shape, -1, np.int32)

        code = np.full(g0.shape, 0, dtype=np.int32)
        if (False & ('p_member' in catalog.data)):
            prob = catalog.data['p_member'][mask]

            top_pri = np.maximum(np.floor((i0 - 16)/(21.5 - 16) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((i0 - 16)/(21.5 - 16) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness
            
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
        else:
            w = (g0 > 18) & (g0 < 23)
            priority[w] = 9

            #young, blue main sequence stars
            w = selection & (g0 > 18) & (g0 < 23) & (gi0 <= 0.0)
            priority[w] = np.minimum(np.maximum(8 - np.rint(4*(g0[w]-23)/(18-23)).astype(int), 0), 8)

            #red giants
            w = selection & (g0 > 18) & (g0 < 23) & (gi0 > 0.0)
            priority[w] = np.minimum(np.maximum(4 - np.rint(4*(g0[w]-23)/(18-23)).astype(int), 0), 8)

        # Blue Horizontal Branch
        #w = (g0 > 19.6) & (g0 < 20.2) & (gi0 > -0.5) & (gi0 < 0.2) & (cli <= 0.5) & (clg < 0.5)
        #priority[w] = 6
        #code[w] = 0
        
        # Very bright stars
        w = (g0 <= 18) & (cli <= 0.5) & (clg <= 0.5)
        priority[w] = 9
        code[w] = 1
        
        # Very faint stars with lowest priority
        w = (g0 >= 23.0) & (cli <= 0.5) & (clg <= 0.5)
        priority[w] = 9
        code[w] = 2

        # Possible extended sources, regardless of magnitude
        w = (cli > 0.5) | (clg > 0.5)
        priority[w] = 9
        code[w] = 3

        # Assign minimum priority to non-members based on Gaia proper motions but within the probability cut
        # These are stars typically a bit bluer than the dSph RGB
        if 'pmra' in catalog.data.columns:
            pmra = catalog.data['pmra'][mask]
            pmra_err = catalog.data['err_pmra'][mask]
            pmdec = catalog.data['pmdec'][mask]
            pmdec_err = catalog.data['err_pmdec'][mask]

            nonmem = (code == 0) & \
                (np.sqrt((pmra - self.pmra.value) ** 2 / (pmra_err ** 2 + self.pmra_err ** 2) +
                            (pmdec - self.pmdec.value) ** 2 / (pmdec_err ** 2 + self.pmdec_err ** 2)) > 3) & \
                (pmra_err >= 0.0) & (pmdec_err >= 0.0) & ~np.isnan(pmra) & ~np.isnan(pmdec)
            priority[nonmem] = 9

        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * ((i0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)
        keep = (g0 < 23) & (priority <= 9) & (code == 0)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]