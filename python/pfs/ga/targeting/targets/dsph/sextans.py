import astropy.units as u
import numpy as np
from datetime import datetime, timedelta, tzinfo

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color
from pfs.ga.common.io import ObservationSerializer

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection, IsochroneSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ...config.pmap import PMapConfig
from ... import Isochrone
from ..ids import *
from .dsphgalaxy import DSphGalaxy

from ...setup_logger import logger

prioritize_DEIMOS = True

class Sextans(DSphGalaxy):
    def __init__(self):
        ID = 'sex'
        name = 'Sextans'

        pos = [ '10h 13m 02.9s', '-01d 36m 53s' ]           # Paturel et al. (2003) Hyperleda
        rad = 52.2 * u.arcmin                               # Forbes et al. (2008)
        DM, DM_err = 19.67, 0.1                             # McConnachie et al. (2012) and references therein
        pm = [ -0.260, 0.100 ] * u.mas / u.yr               # Walker et al. (2008)
        pm_err = [ 0.410, 0.440 ] * u.mas / u.yr            # 
        RV, RV_err = (-224.2, 0.1) * u.kilometer / u.second # McConnachie et al. (2012)

        # ra0 = [ 14.5, 15.1, 15.5, 15.0, 16.4, 15.06, 13.7, 14.9 ] * u.deg
        # dec0 = [ -33.7, -33.4, -33.7, -34.1, -33.9, -33.0, -33.55, -34.5 ] * u.deg
        # pa0 = [ 30, 30, 30, 30, 30, 30, 30, 30 ] * u.deg

        # pointings = {
        #     SubaruPFI: [ Pointing((ra, dec), posang=pa) for ra, dec, pa in zip(ra0, dec0, pa0) ]
        # }

        # Generate the poitings algorithmically
        pointings = {
            SubaruPFI: [
                # Inner pointings along the minor axis
                Pointing.from_relative_pos(pos, sep=0.35, dir=48, posang=220, stage=1, priority=2),
                Pointing.from_relative_pos(pos, sep=-0.35, dir=48, posang=220, stage=1, priority=2),

                # Rotated outer pointings along the minor axis
                Pointing.from_relative_pos(pos, sep=0.85, dir=50, posang=190, stage=3, priority=8),
                Pointing.from_relative_pos(pos, sep=-0.85, dir=50, posang=190, stage=3, priority=8),

                # Inner pointings along the major axis
                Pointing.from_relative_pos(pos, sep=0.45, dir=140, posang=220, stage=0, priority=1),
                Pointing.from_relative_pos(pos, sep=-0.45, dir=140, posang=220, stage=0, priority=1),

                # Outer pointings along the major axis
                Pointing.from_relative_pos(pos, sep=1.05, dir=140, posang=190, stage=2, priority=4),
                Pointing.from_relative_pos(pos, sep=-1.05, dir=140, posang=190, stage=2, priority=4),
            ]
        }

        super().__init__(ID, name, ID_PREFIX_SEXTANS,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        hsc = SubaruHSC.photometry()
        self._hsc_cmd = CMD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i2']]), limits=(-1, 4)),
            MagnitudeAxis(hsc.magnitudes['g'], limits=(15.5, 24.5))
        ])
        self._hsc_ccd = CCD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i2']]), limits=(-1, 4)),
            ColorAxis( Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), limits=(-0.5, 0.5))
        ])

        gaia = Gaia.photometry()
        self._gaia_cmd = CMD([
            ColorAxis(Color([gaia.magnitudes['bp'], gaia.magnitudes['rp']]), limits=(0, 3)),
            MagnitudeAxis(gaia.magnitudes['g'], limits=(11, 22))
        ])

    def get_pmap_config(self):
        config = PMapConfig(
            cut_nb = True,
            keep_blue = True,
            extents = [[0.1, 2.0], [17.0, 23.5]],
            merge_list = [np.s_[:10], np.s_[10:]]
        )

        return config

    def get_filter_map(self):
        """
        Return a dictionary that maps between filter names used internally and the actual
        filter names that are to be written into the exported target lists and design files.

        This is just a final hack because propagating the new filter names through the stack
        can take a significant effort. 
        """

        # Ursa Minor was imaged by Sakurako Okamoto with i band in 2015. That makes it i_hsc
        # instead of i2_hsc.
        filter_map = {
            # 'r_hsc': 'r2_hsc',
            'i_hsc': 'i_old_hsc',
        }

        return filter_map

    def get_text_observation_reader(self, instrument=SubaruHSC):
        if instrument == SubaruHSC:
            mags = ['g', 'i', 'nb515']
            new_mags = {'g': 'g', 'i': 'i2', 'nb515': 'nb515'}
            ext = ['g', 'i', 'nb515']

            reader = ObservationSerializer(format='.csv')
            reader.append_photometry(SubaruHSC.photometry())
            reader.column_map = {
                # 'ID': 'objid',
                'RA': 'RA',
                'Dec': 'Dec',
            }
            reader.column_names = ['RA', 'Dec', 'xi', 'eta']

            for m in mags:
                reader.column_map[f'{m[0]}psf'] = f'obs_hsc_{new_mags[m]}'
                reader.column_map[f'{m[0]}psferr'] = f'err_hsc_{new_mags[m]}'

            for m in ext:
                reader.column_map[f'a_{m[0]}'] = f'ext_hsc_{new_mags[m]}'

            # These have to be done is separate loops because column order matters!

            for m in mags:
                reader.column_names.append(f'{m[0]}psf')

            for m in mags:
                reader.column_names.append(f'{m[0]}psferr')

            for m in mags:
                reader.column_names.append(f'cl{m[0]}')

            reader.column_names.append('EBV')

            for m in ext:
                reader.column_names.append(f'a_{m[0]}')

            for m in ext:
                reader.column_names.append(f'{m}psf0')

            def filter(df):
                ff = None
                for m in mags:
                    if m != 'nb515':
                        f = df[f'cl{m[0]}'] < 0.1
                        if ff is None:
                            ff = f
                        else:
                            ff = ff | f
                return ff

            reader.filter = filter
            reader.kwargs = dict(
                delimiter=',',
                skiprows=1
            )

            return reader
        else:
            raise NotImplementedError()

    def get_nb_selection_mask(self, catalog: Catalog, observed=None, mask=None):
        ccd = self._hsc_ccd

        return (
            ColorSelection(ccd.axes[0], 0.12, 0.5).apply(catalog, observed=observed, mask=mask)

            | ColorSelection(ccd.axes[1], 0.1, None).apply(catalog, observed=observed, mask=mask)
            & ColorSelection(ccd.axes[0], None, 1.65).apply(catalog, observed=observed, mask=mask)
            
            | LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed, mask=mask)
        )

    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self._hsc_cmd
        ccd = self._hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.12, 2.0).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= self.get_nb_selection_mask(catalog, observed=observed)

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

        # Make sure only point sources are selected (applies to observations only)
        if 'clg' in catalog.data and 'cli' in catalog.data:
            mask &= (catalog.data['clg'] < 0.1) & (catalog.data['cli'] < 0.1)

        return mask
    
    def assign_priorities(self, catalog: Catalog, mask=None, isogrid=None):
        """
        Assign priority classes based on photometry
        """

        mask = mask.copy() if mask is not None else np.full(catalog.shape[0], True, dtype=bool)

        logger.info(f'Assigning prioritis to {self.name} stars.')
        logger.info(f'HSC catalog size: {catalog.shape[0]}, using {mask.sum()} unmasked stars.')

        hsc = SubaruHSC.photometry()
        cmd = self._hsc_cmd
        ccd = self._hsc_ccd

        g0, _ = catalog.get_magnitude(hsc.magnitudes['g'], observed=True, dered=True, mask=mask)
        i0, _ = catalog.get_magnitude(hsc.magnitudes['g'], observed=True, dered=True, mask=mask)
        gi0, _ = catalog.get_color(Color([hsc.magnitudes['g'], hsc.magnitudes['i2']]), observed=True, dered=True, mask=mask)
        gn0, _ = catalog.get_color(Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), observed=True, dered=True, mask=mask)

        # Exclude very bright and very faint stars in case they accidentally
        # got into the sample
        keep = mask & (16 <= g0) & (g0 <= 23) & \
                      (16 <= i0) & (i0 < 23)

        # Exclude any extended sources
        clg = catalog.data['clg'][mask]
        cli = catalog.data['cli'][mask]
        keep &= (cli < 0.5) & (clg < 0.5)
        
        p_member = catalog.data['p_member'][mask]
        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10 ** ((i0 - 19.0) / 2.5)) + 1).astype(int), 6), 1)

        # Priorities
        priority = np.full_like(p_member, -1, np.int32)

        # Everything without membership probability
        w9 = np.isnan(p_member) | (p_member == 0.0)
        priority[w9] = 9

        # Priority 0: bright likely members, this will be extended with DEIMOS targets
        w0 = np.isfinite(p_member) & (p_member > 0.8) & (g0 < 22.5)
        priority[w0] = 0
        logger.info(f'{(keep & w0).sum()} {self.name} stars are marked as priority 0')

        # Make cuts based on brightness

        # # Priority 1:
        # w1 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.5) & (g0 < 20.0)
        # priority[w1] = 1
        # logger.info(f'{(keep & w1).sum()} {self.name} stars are marked as priority 1')

        # # Priority 2:
        # w2 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.5) & (g0 < 21.5)
        # priority[w2] = 2
        # logger.info(f'{(keep & w2).sum()} {self.name} stars are marked as priority 2')

        # # Priority 3:
        # w3 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.5)
        # priority[w3] = 3
        # logger.info(f'{(keep & w3).sum()} {self.name} stars are marked as priority 3')

        # Make cuts based on membership

        # Priority 1:
        w1 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.7) & (g0 < 22.5)
        priority[w1] = 1
        logger.info(f'{(keep & w1).sum()} {self.name} stars are marked as priority 1')

        # Priority 2:
        w2 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0) & (g0 < 22.5)
        priority[w2] = 2
        logger.info(f'{(keep & w2).sum()} {self.name} stars are marked as priority 2')

        # Special cuts to include blue stars

        # Priority 3:

        # Blue Horizontal Branch
        wHB = (priority == 9) & (g0 > 19.4) & (g0 < 20.2) & (gi0 > -0.5) & (gi0 < 0.2)
        priority[wHB] = 3
        logger.info(f'{(keep & wHB).sum()} {self.name} BHB stars are marked as priority 3')

        # Potential AGB stars
        # These are stars that are bluer than the RGB but has no membership estimate
        # because we have no reliable models for them
        if isogrid is not None:
            iso_blue = Isochrone()
            iso_blue.from_isogrid(hsc, isogrid, Fe_H=-2.0, log_t=10.111, DM=19.2)
            iso_sel = IsochroneSelection(iso_blue, self._hsc_cmd.axes, selection_axis=0,
                                         selection_direction='-', DM=19.2, error_sigma=[0, 0])
            wAGB = iso_sel.apply(catalog, mask=mask)
            wAGB &= (priority == 9) & (g0 > 17.25) & (g0 < 19.8) & (gi0 > -0.5) & \
                    (gn0 > 0.25 * gi0 - 0.15) & \
                    (g0 > 19.25 - 1.5 * gi0)
            priority[wAGB] = 3

            logger.info(f'{(keep & wAGB).sum()} potential {self.name} AGB stars are marked as priority 3')

        # Stars around the tip of the RGB but red of the halo edge
        wT = (priority == 9) & \
             (16.8 <= g0) & (g0 <= 18.5) & (0.75 <= gi0) & (gi0 <= 1.8) & \
             self.get_nb_selection_mask(catalog, observed=True, mask=mask)
        priority[wT] = 3
        logger.info(f'{(keep & wT).sum()} potential {self.name} bright RGB stars are marked as priority 3')

        # Priority 8: - faint likely members, there's a lot of them
        w8 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0)
        priority[w8] = 8
        logger.info(f'{(keep & w8).sum()} {self.name} stars are marked as priority 8')

        # Priority 4:

        # Blue Stragglers
        x1, x2, x3, x4 = -0.7, -0.2, -0.2, 0.3
        y1, y2, y3, y4 = 20.8, 20.8, 23.0, 23.0
        wBS = (priority == 9) & \
              (g0 > y1) & (g0 < y3) & ((g0 - y1) < ((y3 - y1) / (x3 - x1) * (gi0 - x1))) & \
              ((g0 - y2) > ((y4 - y2) / (x4 - x2) * (gi0 - x2)))
        priority[wBS] = 4
        logger.info(f'{(keep & wBS).sum()} {self.name} Blue Straggler stars are marked as priority 4')

        # Assign minimum priority to non-members based on Gaia proper motions but within the probability cut
        # These are stars typically a bit bluer than the dSph RGB
        if 'pmra' in catalog.data.columns:
            pmra = catalog.data['pmra'][mask]
            pmra_err = catalog.data['err_pmra'][mask]
            pmdec = catalog.data['pmdec'][mask]
            pmdec_err = catalog.data['err_pmdec'][mask]

            nonmem = (priority >= 0) & (p_member > 0) & \
                (np.sqrt((pmra - self.pmra.value) ** 2 / (pmra_err ** 2 + self.pmra_err ** 2) +
                         (pmdec - self.pmdec.value) ** 2 / (pmdec_err ** 2 + self.pmdec_err ** 2)) > 3) & \
                (pmra_err >= 0.0) & (pmdec_err >= 0.0) & ~np.isnan(pmra) & ~np.isnan(pmdec)
            priority[nonmem] = 9

            logger.info(f'{(keep & nonmem).sum()} GAIA stars with high pm are demoted to priority 9')

        # If a target has been previously observed by DEIMOS (Kirby et al. 2010), then set its priority to 0.
        # TODO: move this to the netflow config
        if prioritize_DEIMOS:
            import os
            from astropy.coordinates import SkyCoord
            from astropy.io import fits

            fn = os.path.expandvars('$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/umi_moogify_member.fits.gz')
            hdul = fits.open(fn)
            deimos = hdul[1].data
            c_deimos = SkyCoord(deimos['RA'] * u.degree, deimos['DEC'] * u.degree)
            c_hsc = catalog.get_skycoords()
            idx, d2d, _ = c_deimos.match_to_catalog_sky(c_hsc)
            priority[idx] = 0

            logger.info(f'{idx.size} DEIMOS stars are marked as priority 0')

        # Only keep stars with valid priority
        keep &= (priority >= 0) & (priority <= 9)

        catalog.data['priority'] = -1
        catalog.data.loc[keep, 'priority'] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data.loc[keep, 'exp_time'] = exp_time[keep]