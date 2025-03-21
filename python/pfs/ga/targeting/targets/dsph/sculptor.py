import astropy.units as u
from datetime import datetime, timedelta, tzinfo

from ...util.args import *
from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from ...photometry import Photometry, Magnitude, Color
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ..ids import *
from .dsphgalaxy import DSphGalaxy

class Sculptor(DSphGalaxy):
    def __init__(self):
        ID = 'scl'
        name = 'Sculptor'

        pos = [ '01h 00m 09.4s', '-33d 42m 32s' ]
        # pos = [ 15.0392, -33.7089 ] * u.deg               # Evan
        rad = 120 * u.arcmin
        DM, DM_err = 19.67, 0.1                             # Evan
        # pm, pm_err = [ 0.09, 0.02 ], [ 0.13, 0.13 ]       # Piatek et al. (2006)
        pm = [ 0.101, -0.156 ] * u.mas / u.yr               # Evan
        pm_err = [ 0.003, 0.002 ] * u.mas / u.yr                          
        RV, RV_err = -110.0 * u.kilometer / u.second, 2.0

        ra0 = [ 14.5, 15.1, 15.5, 15.0, 16.4, 15.06, 13.7, 14.9 ] * u.deg
        dec0 = [ -33.7, -33.4, -33.7, -34.1, -33.9, -33.0, -33.55, -34.5 ] * u.deg
        pa0 = [ 30, 30, 30, 30, 30, 30, 30, 30 ] * u.deg

        pointings = {
            SubaruPFI: [ Pointing((ra, dec), posang=pa) for ra, dec, pa in zip(ra0, dec0, pa0) ]
        }

        super().__init__(ID, name, ID_PREFIX_SCULPTOR,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
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

    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self._hsc_cmd
        ccd = self._hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.12, 2.0).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= (
                ColorSelection(ccd.axes[0], 0.12, 0.5).apply(catalog, observed=observed)

                | ColorSelection(ccd.axes[1], 0.1, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], None, 2.0).apply(catalog, observed=observed)
                
                | LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed)
            )

        # Probability-based cut (map) - nonzero membership probability
        if probcut is not None:
            mask &= probcut.apply(catalog, observed=observed)

        # Allow blue
        if blue:
            mask |= (
                ColorSelection(self.ccd.axes[0], None, 0.12).apply(catalog, observed=observed)
            )

        # Always impose faint and bright magnitude cuts
        mask &= MagnitudeSelection(cmd.axes[1], bright, faint).apply(catalog, observed=observed)

        return mask
    
    def assign_priorities(self, catalog: Catalog, mask=None):
        """Assign priority classes based on photometry"""

        # TODO: merge this with Umi
        # TODO: log warnings when decisions are made based on available columns

        mask = mask if mask is not None else np.s_[:]

        g0, i0, gi0 = self._get_hsc_dered_mags_colors(catalog, mask)
        clg = catalog.data['clg'][mask]
        cli = catalog.data['cli'][mask]

        priority = np.full(g0.shape, -1, np.int32)

        if 'p_member' in catalog.data:
            prob = catalog.data['p_member'][mask]
            code = np.full(prob.shape, 0, dtype=np.int32)

            top_pri = np.maximum(np.floor((i0 - 16)/(21.5 - 16) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((i0 - 16)/(21.5 - 16) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness
          
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # Blue Horizontal Branch
            w = (g0 > 19.6) & (g0 < 20.2) & (gi0 > -0.5) & (gi0 < 0.2) & (cli <= 0.5) & (clg < 0.5)
            priority[w] = 6
            code[w] = 0
            
            # Very bright stars, this does nothing because there aren't any of those
            w = (i0 <= 16) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = (i0 >= 21.5) & (cli <= 0.5) & (clg <= 0.5)
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

                nonmem = (code == 0) & (prob > 0) & \
                    (np.sqrt((pmra - self.pmra.value) ** 2 / (pmra_err ** 2 + self.pmra_err ** 2) +
                             (pmdec - self.pmdec.value) ** 2 / (pmdec_err ** 2 + self.pmdec_err ** 2)) > 3) & \
                    (pmra_err >= 0.0) & (pmdec_err >= 0.0) & ~np.isnan(pmra) & ~np.isnan(pmdec)
                priority[nonmem] = 9
        else:
            predicted_y1 = -(0.65 / 0.6) * gi0 + 22.4
            predicted_y2 = -(1.75 / 0.9) * gi0 + 23.5
            predicted_y3 = -(5.75 / 0.3) * gi0 + 33.25
            predicted_y4 = -(5.75 / 0.65) * gi0 + 28.38

            w = (g0 > 16) & (g0 < 23)
            priority[w] = np.rint(11 * (g0[w] - 16) / (23 - 16)).astype(int) + 1
            
            w = (gi0 < 0.6) & ((i0 < predicted_y1) | (i0 > predicted_y2))
            priority[w] = 13
            
            w = (gi0 >= 0.6) & ((i0 < predicted_y3) | (i0 > predicted_y4))
            priority[w] = 13
                
            w = (g0 <= 16)
            priority[w] = 13
            
            w = (g0 >= 23)
            priority[w] = 13
            
            # Cuts on probability of being a point source are already imposed when loading the data file
            # w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            # priority[w] = 14

        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * ((i0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)

        keep = (g0 < 23) & (priority <= 9) & (code == 0)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]