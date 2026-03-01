import numpy as np
import astropy.units as u
from datetime import datetime, timedelta, tzinfo

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...io import ObservationSerializer
from ... import Isochrone
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection, IsochroneSelection
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ...config.pmap import PMapConfig
from ...config.sample import SampleConfig
from ..ids import *
from .dsphgalaxy import DSphGalaxy

from ...setup_logger import logger

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
        
        hsc = SubaruHSC.photometry()
        self._hsc_cmd = CMD([
            ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['r']]), limits=(-0.5, 2)),
            MagnitudeAxis(hsc.magnitudes['g'], limits=(16.5, 24.5))
        ])
        self._hsc_ccd = None
        #self._hsc_ccd = CCD([
        #    ColorAxis(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), limits=(-1, 4)),
        #    ColorAxis( Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), limits=(-0.5, 0.5))
        #])

        gaia = Gaia.photometry()
        self._gaia_cmd = CMD([
            ColorAxis(Color([gaia.magnitudes['bp'], gaia.magnitudes['rp']]), limits=(0, 3)),
            MagnitudeAxis(gaia.magnitudes['g'], limits=(11, 22))
        ])

        cfht = CFHT.photometry()
        self._cfht_cmd = CMD([
            ColorAxis(Color([cfht.magnitudes['g'], cfht.magnitudes['r']]), limits=(-1.5, 3)),
            MagnitudeAxis(cfht.magnitudes['g'], limits=(15.5, 24.5))
        ])

        sdss = SDSS.photometry()
        self._sdss_cmd = CMD([
            ColorAxis(Color([sdss.magnitudes['g'], sdss.magnitudes['r']]), limits=(-1.5, 3)),
            MagnitudeAxis(sdss.magnitudes['g'], limits=(15.5, 24.5))
        ])

    def get_text_observation_reader(self, instrument=SubaruHSC):
        if instrument == SubaruHSC:
            # This is the reader config for the HSC × PS1 × GAIA file from Kohei

            # mags = ['g', 'r']
            # ext = ['g', 'r']

            reader = ObservationSerializer(format='.csv')
            reader.append_photometry(SubaruHSC.photometry())
            reader.append_photometry(PS1.photometry())
            reader.append_photometry(Gaia.photometry())


            reader.column_map = {
                'hsc_targetid': 'objid',
                'hsc_ra': 'RA',
                'hsc_dec': 'Dec',
                
                'hsc_g': 'obs_hsc_g',
                'hsc_r': 'obs_hsc_r',
                'hsc_a_g': 'ext_hsc_g',
                'hsc_a_r': 'ext_hsc_r',

                'ps1_g': 'obs_ps1_g',
                'ps1_r': 'obs_ps1_r',
                'ps1_i': 'obs_ps1_i',
                'ps1_z': 'obs_ps1_z',
                'ps1_y': 'obs_ps1_y',

                # 'gaia_g': 'obs_gaia_g',
                # 'gaia_bp': 'obs_gaia_bp',
                # 'gaia_rp': 'obs_gaia_rp',
            }
            
            #reader.column_names = ['ID', 'RA', 'Dec', 'X', 'Y']

            reader.kwargs = dict(
                delimiter=','
            )

            return reader
        elif instrument == CFHT:
            return CFHT.text_observation_reader(
                mags=['g', 'r'], ext=['g', 'r'])
        else:
            raise NotImplementedError()

    def get_cmd(self, index=0, instrument=SubaruHSC):
        """
        Return the definition of the color-magnitude diagram for the given instrument.
        """

        if instrument == SubaruHSC:
            if index is not None and isinstance(self._hsc_cmd, Iterable):
                cmd = self._hsc_cmd[index]
            else:
                cmd = self._hsc_cmd
        elif instrument == Gaia:
            cmd = self._gaia_cmd
        elif instrument == CFHT:
            cmd = self._sdss_cmd        # Munoz
        else:
            raise NotImplementedError()
            
        return cmd

    def get_pmap_config(self):
        config = PMapConfig(
            cut_nb = False,
            keep_blue = False,
            extents = [[0.0, 1.5], [16.8, 23]],
            merge_list = [np.s_[:10], np.s_[10:]]
        )

        return config

    def get_sample_config(self):
        config = SampleConfig()
        return config
    
    def get_selection_mask(self, catalog: Catalog, nb=False, blue=False, probcut=None, observed=None, bright=17, faint=23):
        """Return true for objects within sharp magnitude cuts."""

        cmd = self._hsc_cmd
        # ccd = self._hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.1, 1.5).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            raise NotImplementedError()

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
    
    # def _get_cfht_dered_mags_colors(self, catalog, mask):
    #     cfht = CFHT.photometry()
    #     [ (g0, _), (r0, _), (gr0, _) ] = catalog.get_diagram_values([
    #                 cfht.magnitudes['g'],
    #                 cfht.magnitudes['r'],
    #                 Color([cfht.magnitudes['g'], cfht.magnitudes['r']])
    #             ], observed=True, mask=mask)
        
    #     raise NotImplementedError()

    #     return g0, r0, gr0

    # def assign_priorities_CFHT(self, catalog: Catalog, mask=None, isogrid=None, isochrones_name_mappings=None):
    #     """Assign priority classes based on photometry"""

    #     mask = mask.copy() if mask is not None else np.full(catalog.shape[0], True, dtype=bool)

    #     logger.info(f'Assigning prioritis to {self.name} stars.')
    #     logger.info(f'HSC catalog size: {catalog.shape[0]}, using {mask.sum()} unmasked stars.')

    #     # CFHT data in Munoz et al. (2018) is given in SDSS magnitudes
    #     cfht = CFHT.photometry()

    #     g0, _ = catalog.get_magnitude(cfht.magnitudes['g'], observed=True, dered=True, mask=mask)
    #     r0, _ = catalog.get_magnitude(cfht.magnitudes['r'], observed=True, dered=True, mask=mask)
    #     gr0, _ = catalog.get_color(Color([cfht.magnitudes['g'], cfht.magnitudes['r']]), observed=True, dered=True, mask=mask)

    #     # Exclude very bright and very faint stars in case they accidentally
    #     # got into the sample
    #     keep = mask & (16 <= g0) & (g0 <= 23) & \
    #                   (16 <= r0) & (r0 < 23) & \
    #                   (-0.4 <= gr0) & (gr0 < 2)
        
    #     p_member = catalog.data['p_member'][mask]
    #     exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10 ** ((r0 - 19.0) / 2.5)) + 1).astype(int), 6), 1)

    #     # We have no narrow-band data for Bootes, so we will not be able to select stars based on that.
    #     # This means that membership probability declines significantly toward the tip of the RGB, hence
    #     # we assign priorities based on brightness rather than membership. None-zero membership is used to
    #     # select the region where RGB stars are likely.

    #     # Priorities
    #     priority = np.full_like(p_member, -1, np.int32)

    #     # Everything without membership probability
    #     w9 = (~np.isfinite(p_member) | (p_member == 0)) & (g0 > 17.5)
    #     priority[w9] = 9

    #     # Priority 0: bright likely members
    #     w0 = np.isfinite(p_member) & (p_member > 0) & (g0 < 20.0)
    #     priority[w0] = 0
    #     logger.info(f'{(keep & w0).sum()} {self.name} stars are marked as priority 0')

    #     # Make cuts based on membership

    #     # Priority 1:
    #     w1 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.5) & (g0 < 22.5)
    #     priority[w1] = 1
    #     logger.info(f'{(keep & w1).sum()} {self.name} stars are marked as priority 1')

    #     # Priority 2:
    #     w2 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0) & (g0 < 22.5)
    #     priority[w2] = 2
    #     logger.info(f'{(keep & w2).sum()} {self.name} stars are marked as priority 2')

    #     # Special cuts to include blue stars

    #     # Priority 3:

    #     # # Blue Horizontal Branch
    #     # wHB = (priority == 9) & (g0 > 19) & \
    #     #       (g0 < 20) & (gr0 > -0.5) & (gr0 < 0.2)
    #     # priority[wHB] = 3
    #     # logger.info(f'{(keep & wHB).sum()} {self.name} BHB stars are marked as priority 3')

    #     # # Potential AGB stars
    #     # # These are stars that are bluer than the RGB but has no membership estimate
    #     # # because we have no reliable models for them
    #     # if isogrid is not None:
    #     #     iso_blue = Isochrone()
    #     #     iso_blue.from_isogrid(cfht, isogrid, Fe_H=-2.5, log_t=10.1, DM=19.11)
    #     #     iso_sel = IsochroneSelection(iso_blue, self._cfht_cmd.axes, selection_axis=0,
    #     #                                  selection_direction='-', DM=19.11, error_sigma=[0, 0])
    #     #     wAGB = iso_sel.apply(catalog, mask=mask)
    #     #     wAGB &= (priority == 9) & (g0 > 17.25) & (g0 < 20.6) & (gr0 > -0.5) & \
    #     #             (g0 > 17)
    #     #             # (g0 > 19.15 - 1.5 * gr0)
    #     #             # (gn0 > 0.25 * gr0 - 0.15) & \
    #     #     priority[wAGB] = 3

    #     #     logger.info(f'{(keep & wAGB).sum()} potential {self.name} AGB stars are marked as priority 3')

    #     # Blue stars bluer than the RGB
    #     # This should include BHB stars and blue stragglers
    #     if isogrid is not None:
    #         iso_blue = Isochrone()
    #         iso_blue.from_isogrid(cfht, isogrid, Fe_H=-2.5, log_t=10.1, DM=19.11)
    #         iso_sel = IsochroneSelection(iso_blue, self._cfht_cmd.axes, selection_axis=0,
    #                                      selection_direction='-', DM=19.11, error_sigma=[0, 0])
    #         wAGB = iso_sel.apply(catalog, mask=mask)
    #         wAGB &= (priority == 9) & (g0 > 17.25) & (g0 < 22) & (gr0 > -0.5)
    #         priority[wAGB] = 3

    #         logger.info(f'{(keep & wAGB).sum()} potential {self.name} AGB stars are marked as priority 3')

    #     # Stars around the tip of the RGB but red of the halo edge
    #     # wT = (priority == 9) & \
    #     #      (16.8 <= g0) & (g0 <= 18.5) & (0.75 <= gr0) & (gr0 <= 1.8) & \
    #     #      self.get_nb_selection_mask(catalog, observed=True, mask=mask)
    #     # priority[wT] = 3
    #     # logger.info(f'{(keep & wT).sum()} potential {self.name} bright RGB stars are marked as priority 3')

    #     # Priority 8: - faint likely members, there's a lot of them
    #     w8 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0)
    #     priority[w8] = 8
    #     logger.info(f'{(keep & w8).sum()} {self.name} stars are marked as priority 8')

    #     # Priority 4:

    #     # Blue Stragglers
    #     x1, x2, x3, x4 = -1.0, -0.1, -0.6, 0.2
    #     y1, y2, y3, y4 = 20.8, 20.8, 23.0, 23.0
    #     wBS = (priority == 9) & \
    #           (g0 > y1) & (g0 < y3) & ((g0 - y1) < ((y3 - y1) / (x3 - x1) * (gr0 - x1))) & \
    #           ((g0 - y2) > ((y4 - y2) / (x4 - x2) * (gr0 - x2)))
    #     priority[wBS] = 4
    #     logger.info(f'{(keep & wBS).sum()} {self.name} Blue Straggler stars are marked as priority 4')

    #     # Assign minimum priority to non-members based on Gaia proper motions but within the probability cut
    #     # These are stars typically a bit bluer than the dSph RGB
    #     if 'pmra' in catalog.data.columns:
    #         pmra = catalog.data['pmra'][mask]
    #         pmra_err = catalog.data['err_pmra'][mask]
    #         pmdec = catalog.data['pmdec'][mask]
    #         pmdec_err = catalog.data['err_pmdec'][mask]

    #         nonmem = (priority >= 0) & (p_member > 0) & \
    #             (np.sqrt((pmra - self.pmra.value) ** 2 / (pmra_err ** 2 + self.pmra_err ** 2) +
    #                      (pmdec - self.pmdec.value) ** 2 / (pmdec_err ** 2 + self.pmdec_err ** 2)) > 3) & \
    #             (pmra_err >= 0.0) & (pmdec_err >= 0.0) & ~np.isnan(pmra) & ~np.isnan(pmdec)
    #         priority[nonmem] = 9

    #         logger.info(f'{(keep & nonmem).sum()} GAIA stars with high pm are demoted to priority 9')

    #     # Only keep stars with valid priority
    #     keep &= (priority >= 0) & (priority <= 9)

    #     catalog.data['priority'] = -1
    #     catalog.data.loc[keep, 'priority'] = priority[keep]

    #     catalog.data['exp_time'] = np.nan
    #     catalog.data.loc[keep, 'exp_time'] = exp_time[keep]

    def assign_priorities(self, catalog: Catalog, mask=None, isogrid=None, isochrones_name_mappings=None):
        """Assign priority classes based on photometry"""

        mask = mask.copy() if mask is not None else np.full(catalog.shape[0], True, dtype=bool)

        logger.info(f'Assigning prioritis to {self.name} stars.')
        logger.info(f'HSC catalog size: {catalog.shape[0]}, using {mask.sum()} unmasked stars.')

        hsc = SubaruHSC.photometry()

        g0, _ = catalog.get_magnitude(hsc.magnitudes['g'], observed=True, dered=True, mask=mask)
        r0, _ = catalog.get_magnitude(hsc.magnitudes['r'], observed=True, dered=True, mask=mask)
        gr0, _ = catalog.get_color(Color([hsc.magnitudes['g'], hsc.magnitudes['r']]), observed=True, dered=True, mask=mask)

        # Exclude very bright and very faint stars in case they accidentally
        # got into the sample
        keep = mask & (16 <= g0) & (g0 <= 23) & \
                      (16 <= r0) & (r0 < 23) & \
                      (-0.4 <= gr0) & (gr0 < 2)
        
        p_member = catalog.data['p_member'][mask]
        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10 ** ((r0 - 19.0) / 2.5)) + 1).astype(int), 6), 1)

        # We have no narrow-band data for Bootes, so we will not be able to select stars based on that.
        # This means that membership probability declines significantly toward the tip of the RGB, hence
        # we assign priorities based on brightness rather than membership. None-zero membership is used to
        # select the region where RGB stars are likely.

        # Priorities
        priority = np.full_like(p_member, -1, np.int32)

        # Everything without membership probability
        w9 = (~np.isfinite(p_member) | (p_member == 0)) & (g0 > 17.5)
        priority[w9] = 9

        # Priority 0: bright likely members
        w0 = np.isfinite(p_member) & (p_member > 0) & (g0 < 20.0)
        priority[w0] = 0
        logger.info(f'{(keep & w0).sum()} {self.name} stars are marked as priority 0')

        # Make cuts based on membership

        # Priority 1:
        w1 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.5) & (g0 < 22.5)
        priority[w1] = 1
        logger.info(f'{(keep & w1).sum()} {self.name} stars are marked as priority 1')

        # Priority 4:
        w2 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0) & (g0 < 22.5)
        priority[w2] = 4
        logger.info(f'{(keep & w2).sum()} {self.name} stars are marked as priority 2')

        # Special cuts to include blue stars

        # Priority 2:

        # Blue Horizontal Branch
        wHB = (priority == 9) & \
              (g0 > 19.2) & (g0 < 20) & \
              (gr0 > -0.5) & (gr0 < 0.1)
        priority[wHB] = 2
        logger.info(f'{(keep & wHB).sum()} {self.name} BHB stars are marked as priority 3')

        # Priority 3:
        # Blue Stragglers
        wBS = (priority == 9) & \
              (21 <= g0) & (g0 <= 23) & \
              (-0.5 <= gr0) & (gr0 <= 0.2)
        priority[wBS] = 3
        logger.info(f'{(keep & wBS).sum()} {self.name} Blue Straggler stars are marked as priority 4')

        # Priority 5: - faint likely members, there's a lot of them
        w8 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.0)
        priority[w8] = 5
        logger.info(f'{(keep & w8).sum()} {self.name} stars are marked as priority 8')

        # Priority 6:

        # Blue stars bluer than the RGB
        if isogrid is not None:
            iso_blue = Isochrone()
            iso_blue.from_isogrid(hsc, isogrid, Fe_H=-2.5, log_t=10.1, DM=19.11)
            iso_sel = IsochroneSelection(iso_blue, self._hsc_cmd.axes, selection_axis=0,
                                         selection_direction='-', DM=19.11, error_sigma=[0, 0])
            wAGB = iso_sel.apply(catalog, mask=mask)
            wAGB &= (priority == 9) & (g0 > 17.2) & (g0 < 22) & (gr0 > -0.5)
            priority[wAGB] = 6

            logger.info(f'{(keep & wAGB).sum()} potential {self.name} AGB stars are marked as priority 3')        

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

        # Only keep stars with valid priority
        keep &= (priority >= 0) & (priority <= 9)

        catalog.data['priority'] = -1
        catalog.data.loc[keep, 'priority'] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data.loc[keep, 'exp_time'] = exp_time[keep]