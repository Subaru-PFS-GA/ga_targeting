import numpy as np
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
from .dsphgalaxy import DSphGalaxy

class Fornax(DSphGalaxy):
    def __init__(self):
        ID = 'fornax'
        name = 'Fornax'
        pos = [ 39.9971, -34.4492 ] * u.deg                     # Evan
        rad = 240 * u.arcmin
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

        super().__init__(ID, name,
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

    def get_netflow_config(self):
        config = NetflowConfig.default()

        config.field = FieldConfig(
            key = self.ID,
            name = self.name,
            center = self.get_center(),
            arms = 'bmn',
            nvisits = 1,
            exp_time = 6 * 30 * 60.,        # 3 hr total
            obs_time = datetime(2024, 11, 21, 0, 0, 0) + timedelta(hours=10),
            resolution = 'm',
        )

        config.pointings = [ PointingConfig(p.ra, p.dec, p.posang) for p in self.get_pointings(SubaruPFI) ]

        return config
        
    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.0):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self._hsc_cmd
        ccd = self._hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.12, 2.5).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= (
                ColorSelection(ccd.axes[0], None, 0.0).apply(catalog, observed=observed)

                | LinearSelection(ccd.axes, [0.25, 1.0], 0.2, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], 0.1, 1.5).apply(catalog, observed=observed)
                
                | LinearSelection(ccd.axes, [-0.33, 1.0], -0.66, None).apply(catalog, observed=observed)
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
    
    def assign_priorities(self, catalog: Catalog, mask=None):
        """Assign priority classes based on photometry"""

        mask = mask if mask is not None else np.s_[:]

        g0, i0, gi0 = self._get_hsc_dered_mags_colors(catalog, mask)
        clg = catalog.data['clg'][mask]
        cli = catalog.data['cli'][mask]

        priority = np.full(g0.shape, -1, np.int32)

        if 'p_member' in catalog.data:
            prob = catalog.data['p_member'][mask]
            code = np.full(prob.shape, 0, dtype=np.int32)

            top_pri = np.maximum(np.floor((g0 - 16)/(23 - 16) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((g0 - 16)/(23 - 16) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness
          
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # blue stars
            w = (gi0 < 0.12) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 0
            
            # Blue Horizontal Branch
            w = (i0 > 20.8) & (i0 < 21.6) & (gi0 > -0.5) & (gi0 < 0.4) & (cli <= 0.5) & (clg < 0.5)
            priority[w] = 6
            code[w] = 0
            
            # Very bright stars, this does nothing because there aren't any of those
            w = (i0 <= 16) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = (g0 >= 23) & (cli <= 0.5) & (clg <= 0.5)
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
