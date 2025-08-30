import numpy as np
import astropy.units as u
from datetime import datetime, timedelta, tzinfo

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ..ids import *
from .dsphgalaxy import DSphGalaxy

class Bootes(DSphGalaxy):
    def __init__(self):
        ID = 'bootes'
        name = 'Bootes I'

        pos = [ 210.025, 14.5 ] * u.deg                       # Evan
        rad = 120 * u.arcmin
        DM, DM_err = 19.11, 0.008                             # Oakes et al. (2022)
        pm = [ -0.387, -1.064 ] * u.mas / u.yr                # Evan
        pm_err = [ 0.122, 0.098 ] * u.mas / u.yr
        RV, RV_err = (99.0, 2.1) * u.kilometer / u.second     # Simbad

        # SSP pointings
        ra0 = [ 210.5, 209.6, 210.1, 210.1 ] * u.deg
        dec0 = [ 14.5, 14.5, 14.15, 14.8 ] * u.deg
        pa0 = [ 30, 30, 30, 30 ] * u.deg

        # Engineering run
        # ra0 = [ 210.025 ] * u.deg
        # dec0 = [ 14.5 ] * u.deg
        # pa0 = [ 30 ] * u.deg

        pointings = {
            SubaruPFI: [ Pointing((ra, dec), posang=pa) for ra, dec, pa in zip(ra0, dec0, pa0) ]
        }

        super().__init__(ID, name, ID_PREFIX_BOOTES,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        hsc = SubaruHSC.photometry()
        self._hsc_cmd = CMD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['r']]), limits=(-1, 4)),
            MagnitudeAxis(hsc.magnitudes['g'], limits=(15.5, 24.5))
        ])
        #self._hsc_ccd = CCD([
        #    ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
        #    ColorAxis( Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), limits=(-0.5, 0.5))
        #])

        gaia = Gaia.photometry()
        self._gaia_cmd = CMD([
            ColorAxis(Color([gaia.magnitudes['bp'], gaia.magnitudes['rp']]), limits=(0, 3)),
            MagnitudeAxis(gaia.magnitudes['g'], limits=(11, 22))
        ])
    
    def get_selection_mask(self, catalog: Catalog, nb=False, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self._hsc_cmd
        #ccd = self._hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.08, 2.0).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            raise NotImplementedError

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
    
    def _get_hsc_dered_mags_colors(self, catalog, mask):
        hsc = SubaruHSC.photometry()
        [ (g0, _), (r0, _), (gr0, _) ] = catalog.get_diagram_values([
                    hsc.magnitudes['g'],
                    hsc.magnitudes['r'],
                    Color([hsc.magnitudes['g'], hsc.magnitudes['r']])
                ], observed=True, mask=mask)
        
        return g0, r0, gr0

    def assign_priorities(self, catalog: Catalog, mask=None, isogrid=None):
        """Assign priority classes based on photometry"""

        mask = mask if mask is not None else np.s_[:]

        g0, r0, gr0 = self._get_hsc_dered_mags_colors(catalog, mask)
        clg = catalog.data['clg'][mask]
        clr = catalog.data['clr'][mask]

        priority = np.full(g0.shape, -1, np.int32)

        if 'p_member' in catalog.data:
            prob = catalog.data['p_member'][mask]
            code = np.full(prob.shape, 0, dtype=np.int32)

            top_pri = np.maximum(np.floor((r0 - 16)/(21.5 - 16) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((r0 - 16)/(21.5 - 16) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness
          
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # Blue Horizontal Branch
            w = (g0 > 19.1) & (g0 < 19.7) & (gr0 > -0.5) & (gr0 < 0.2) & (clr <= 0.5) & (clg < 0.5)
            priority[w] = 6
            code[w] = 0
            
            # Very bright stars, this does nothing because there aren't any of those
            w = (r0 <= 16) & (clr <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = (r0 >= 23) & (clr <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 2

            # Possible extended sources, regardless of magnitude
            w = (clr > 0.5) | (clg > 0.5)
            priority[w] = 9
            code[w] = 3

            # Assign minimum priority to non-members based on Gaia proper motions but within the probability cut
            # These are stars typically a bit bluer than the dSph RGB
            if 'pmra' in catalog.data.columns:
                pmra = catalog.data['pmra'][mask]
                pmra_err = catalog.data['err_pmra'][mask]
                pmdec = catalog.data['pmdec'][mask]
                pmdec_err = catalog.data['err_pmdec'][mask]

                nonmem = (code == 0) & (prob > 0) & \
                    (np.sqrt((pmra - self.pmra.value) ** 2 / (pmra_err ** 2 + self.pmra_err ** 2) +
                             (pmdec - self.pmdec.value) ** 2 / (pmdec_err ** 2 + self.pmdec_err ** 2)) > 3) & \
                    (pmra_err >= 0.0) & (pmdec_err >= 0.0) & ~np.isnan(pmra) & ~np.isnan(pmdec)
                priority[nonmem] = 9
        else:
           raise NotImplementedError

        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)

        keep = (g0 < 23) & (priority <= 9) & (code == 0)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]