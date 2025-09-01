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

class M31(Galaxy):

    FIELDS = [
        # group_name, field_name, ra, dec, posang, stage, priority
        ('sector_0', 'PFS_2',   '00:50:39.13',  '+41:32:16.9', 50.0, 0, 0),
        ('sector_0', 'PFS_3',   '00:50:01.64',  '+40:24:31.2', 50.0, 0, 0),
        ('sector_1', 'PFS_4',   '00:44:27.95',  '+39:41:48.1', 50.0, 0, 0),
        ('sector_2', 'PFS_6',   '00:56:29.68',  '+41:54:14.8', 50.0, 0, 0),
        ('sector_2', 'PFS_8',   '00:54:01.84',  '+39:39:08.5', 50.0, 0, 0),
        ('sector_1', 'PFS_10',   '00:42:49.51',  '+38:53:14.9', 50.0, 0, 0),
        ('sector_0', 'PFS_1',   '00:51:50.50',  '+42:39:59.8', 50.0, 0, 0),
        ('sector_2', 'PFS_7',   '00:55:14.42',  '+40:46:43.3', 50.0, 0, 0),
        ('sector_2', 'PFS_9',   '00:48:23.81',  '+39:16:42.8', 50.0, 0, 0),
        ('sector_1', 'PFS_5',   '00:38:59.79',  '+39:14:15.7', 50.0, 0, 0),
        ('sector_2', 'PFS_13',   '00:58:27.46',  '+38:53:07.1', 50.0, 0, 0),
        ('sector_2', 'PFS_11',   '01:01:01.75',  '+41:04:59.3', 50.0, 0, 0),
        ('sector_2', 'PFS_12',   '00:59:49.28',  '+40:03:07.7', 50.0, 0, 0),
        ('sector_2', 'PFS_14',   '00:52:39.04',  '+38:30:07.3', 50.0, 0, 0),
        ('sector_2', 'PFS_15',   '00:47:33.81',  '+38:09:19.1', 50.0, 0, 0),
    ]
    """
        ('sector_0', 'PFS_1',   '00:34:54'  , '+42:22:57', 0.0, 0, 0),
        ('sector_1', 'PFS_2',   '00:28:56'  , '+43:11:01', 0.0, 0, 0),
        ('sector_1', 'PFS_3',   '00:22:49'  , '+43:57:54', 0.0, 0, 0),
        ('sector_1', 'PFS_6',   '00:41:43'  , '+42:54:25', 0.0, 0, 0),
        ('sector_1', 'PFS_7',   '00:35:47'  , '+43:43:47', 0.0, 0, 0),
        ('sector_1', 'PFS_8',   '00:29:42'  , '+44:32:01', 0.0, 0, 0),
        ('sector_1', 'PFS_9',   '00:23:26'  , '+45:19:01', 0.0, 0, 0),
        ('sector_1', 'PFS_11',  '00:42:45'  , '+44:15:02', 0.0, 0, 0),
        ('sector_1', 'PFS_12',  '00:36:42'  , '+45:04:35', 0.0, 0, 0),
        ('sector_1', 'PFS_13',  '00:37:40'  , '+46:25:20', 0.0, 0, 0),
        ('sector_1', 'PFS_14',  '00:34:04'  , '+41:02:05', 0.0, 0, 0),
        ('sector_1', 'PFS_15',  '00:28:13'  , '+41:50:00', 0.0, 0, 0),
        ('sector_1', 'PFS_16',  '00:22:13'  , '+42:36:45', 0.0, 0, 0),
        ('sector_1', 'PFS_19',  '00:50:18'  , '+40:07:25', 0.0, 0, 0),
        ('sector_1', 'PFS_20',  '00:55:44'  , '+39:15:27', 0.0, 0, 0),
        ('sector_1', 'PFS_21',  '01:01:02'  , '+38:22:34', 0.0, 0, 0),
        ('sector_2', 'PFS_22',  '00:51:28'  , '+41:27:45', 0.0, 0, 0),
        ('sector_1', 'PFS_23',  '00:56:59'  , '+40:35:34', 0.0, 0, 0),
        ('sector_1', 'PFS_24',  '01:02:22'  , '+39:42:28', 0.0, 0, 0),
        ('sector_1', 'PFS_25',  '00:52:41'  , '+42:48:01', 0.0, 0, 0),
        ('sector_1', 'PFS_26',  '00:58:18'  , '+41:55:37', 0.0, 0, 0),
        ('sector_1', 'PFS_27',  '01:03:46'  , '+41:02:17', 0.0, 0, 0),
        ('sector_1', 'PFS_28',  '00:43:43'  , '+39:37:51', 0.0, 0, 0),
        ('sector_1', 'PFS_29',  '00:49:12'  , '+38:47:03', 0.0, 0, 0),
        ('sector_1', 'PFS_30',  '00:54:32'  , '+37:55:17', 0.0, 0, 0),
        ('sector_1', 'PFS_31',  '00:37:13'  , '+39:06:52', 0.0, 0, 0),
        ('sector_1', 'PFS_32',  '00:42:44'  , '+38:17:16', 0.0, 0, 0),
        ('sector_1', 'PFS_33',  '00:48:06'  , '+37:26:39', 0.0, 0, 0),
        ('sector_1', 'M31_003', '00:16:31.7', '+44:43:30', 0.0, 0, 0),
        ('sector_1', 'M31_004', '00:10:04.7', '+45:27:47', 0.0, 0, 0),
        ('sector_1', 'M31_009', '00:10:24.5', '+46:49:07', 0.0, 0, 0),
        ('sector_1', 'M31_022', '00:16:04.2', '+43:22:15', 0.0, 0, 0),
        ('sector_1', 'M31_023', '00:09:45.8', '+44:06:28', 0.0, 0, 0),
        ('sector_1', 'PFS_34',  '00:06:10'  , '+47:35:07', 0.0, 0, 0),
        ('sector_1', 'PFS_35',  '00:52:10'  , '+44:25:07', 0.0, 0, 0),
        ('sector_1', 'PFS_36',  '00:48:10'  , '+43:30:07', 0.0, 0, 0),
        ('sector_1', 'PFS_37',  '00:47:06'  , '+42:16:07', 0.0, 0, 0),
        ('sector_1', 'PFS_38',  '00:47:05'  , '+40:40:07', 0.0, 0, 0),
        ('sector_1', 'PFS_39',  '00:39:30'  , '+41:50:07', 0.0, 0, 0),
        ('sector_1', 'PFS_40',  '00:37:40'  , '+40:05:04', 0.0, 0, 0),
        ('sector_1', 'PFS_41',  '00:32:40'  , '+40:05:04', 0.0, 0, 0),
        ('sector_1', 'PFS_42',  '00:31:10'  , '+39:05:04', 0.0, 0, 0),
    """

    def __init__(self, sector=None, field=None):
    
        ID = 'm31'
        name = 'M31'

        pos = [ '00h 42m 44.33s', '+41d 16m 07.5s' ]
        rad = 6 * u.deg
        DM, DM_err = 24.407, 0.032   #Siyang Li et al. (2021, ApJ, 920, 84)
        pm = [ -0.0533, -0.0104 ] * u.mas / u.yr    #Sohn et al. (2012)
        pm_err = [ 0.0246, 0.0244 ] * u.mas / u.yr
        RV, RV_err = (-287.2, 8.0) * u.kilometer / u.second

        pointings_by_sector = defaultdict(dict)
        for sc, fl, ra, dec, pa0, st, pri in self.FIELDS:
            pointings_by_sector[sc][fl] = Pointing(
                Angle(ra, unit=u.hourangle),
                Angle(dec, unit=u.deg),
                posang=pa0, priority=pri, stage=st,
                exp_time=30*60, nvisits=10)
        
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

        super().__init__(ID, name, ID_PREFIX_M31,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        # CMD and CCD definitions with color and magnitude limits for M31

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

    def get_text_observation_reader(self, instrument=SubaruHSC):
        if instrument == SubaruHSC:
            return SubaruHSC.text_observation_reader_m31(
                mags=['g', 'i', 'nb515'],
                ext=['g', 'i', 'nb515'])
        else:
            raise NotImplementedError()

    def get_photometry(self, instrument=SubaruHSC):
        """
        Return the photometric system for the given instrument.
        """
        return instrument.photometry()

    def get_cmd(self, instrument=SubaruHSC):
        """
        Return the definition of the color-magnitude diagram for the given instrument.
        """
        if instrument == SubaruHSC:
            cmd = self.__hsc_cmd
        elif instrument == Gaia:
            cmd = self.__gaia_cmd
        else:
            raise NotImplementedError()
            
        return cmd
    
    def get_ccd(self, instrument=SubaruHSC):
        """
        Return the definition of the color-color diagram for the given instrument.
        """
        if instrument == SubaruHSC:
            ccd = self.__hsc_ccd
        else:
            raise NotImplementedError()

        return ccd

    def get_pmap_config(self):
        config = PMapConfig(
            cut_nb = True,
            keep_blue = False,
            extents = [[0.1, 3.0], [17.0, 24.5]],
            merge_list = [np.s_[:10], np.s_[10:]]
        )

        return config

    def get_sample_config(self):
        config = SampleConfig()
        return config
    
    def get_netflow_config(self):
        config = NetflowConfig.default()

        config.field = self.get_field_config()
        config.pointings = [ PointingConfig.from_pointing(p) for p in self.get_pointings(SubaruPFI) ]

        return config

    def get_field_config(self):
        """
        Return the default field configuration for M31.
        """

        # TODO: it doesn't check visibility at abs_time

        return FieldConfig(
            key = self.ID,
            name = self.name,
            id_prefix = self.id_prefix,
            center = PointingConfig.from_pointing(self.get_center()),
            arms = 'bmn',
            nvisits = 10,
            exp_time = 30 * 60.,        # 5 hr total
            obs_time = datetime(2025, 9, 13, 23, 0, 0) + timedelta(hours=10),
            resolution = 'm',
        )

    def get_filter_map(self):
        """
        Return a dictionary that maps between filter names used internally and the actual
        filter names that are to be written into the exported target lists and design files.

        This is just a final hack because propagating the new filter names through the stack
        can take a significant effort. 
        """

        # What HSC filters were used to observe M31?
        filter_map = {
            # 'r_hsc': 'r2_hsc',
            # 'i_hsc': 'i2_hsc',
        }

        return filter_map
    
    def lookup_ml_prob(self, catalog: Catalog, mask=None):
        """
        Returns the color-diagram-based membership probability based on
        machine learning-based probability (K. Ding et al. 2025).
        """

        ml_prob = catalog.data['ml_prob']
        mask_member = np.full_like(ml_prob, True, dtype=bool)
        mask_member[np.isnan(ml_prob)] = False

        return ml_prob, mask_member

    def assign_probabilities(self, catalog, pmap, population_id=-1, mask=None):
         # Membership probability
        ml_member, ml_member_mask = self.lookup_ml_prob(catalog, mask=mask)

        catalog.data['p_member'] = np.nan
        catalog.data.loc[ml_member_mask,'p_member'] = ml_member[ml_member_mask]

    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=19.0, faint=30.0):
        """Return true for objects within sharp magnitude cuts."""

        cmd = self.__hsc_cmd
        ccd = self.__hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.95, 4.0).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= (
                ColorSelection(ccd.axes[0], 0.12, 0.5).apply(catalog, observed=observed)

                | ColorSelection(ccd.axes[1], 0.05, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], None, 4.0).apply(catalog, observed=observed)
                
                #| LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed)
                # & LinearSelection(ccd.axes, [1.0, 7.0], 1.8, None).apply(catalog, observed=observed)
            )

        # Probability-based cut (map) - nonzero membership probability
        #if probcut is not None:
        #    mask &= probcut.apply(catalog, observed=observed)

        # Allow blue
        if blue:
            mask |= (
                ColorSelection(ccd.axes[0], None, 0.12).apply(catalog, observed=observed)
            )

        # Always impose faint and bright magnitude cuts
        mask &= MagnitudeSelection(cmd.axes[1], bright, faint).apply(catalog, observed=observed)

        return mask
    
    def assign_priorities(self, catalog: Catalog, mask=None, isogrid=None):
        """
        Assign priority classes based on photometry
        """

        mask = mask.copy() if mask is not None else np.full(catalog.shape[0], True, dtype=bool)

        logger.info(f'Assigning priorities to {self.name} stars.')
        logger.info(f'HSC catalog size: {catalog.shape[0]}, using {mask.sum()} unmasked stars.')

        hsc = SubaruHSC.photometry()
        # cmd = self._hsc_cmd
        # ccd = self._hsc_ccd

        g0, _ = catalog.get_magnitude(hsc.magnitudes['g'], observed=True, dered=True, mask=mask)
        i0, _ = catalog.get_magnitude(hsc.magnitudes['i'], observed=True, dered=True, mask=mask)
        # gi0, _ = catalog.get_color(Color([hsc.magnitudes['g'], hsc.magnitudes['i']]), observed=True, dered=True, mask=mask)
        # gn0, _ = catalog.get_color(Color([hsc.magnitudes['g'], hsc.magnitudes['nb515']]), observed=True, dered=True, mask=mask)

        # Exclude very bright and very faint stars in case they accidentally
        # got into the sample
        keep = mask & (16 <= g0) & (g0 <= 30.0) & \
                      (16 <= i0) & (i0 < 30.0)

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

        # Priority 1:
        w1 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.95) & (i0 < 22.25)
        priority[w1] = 1
        logger.info(f'{(keep & w1).sum()} {self.name} stars are marked as priority 1')

        # Priority 2:
        w2 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.95) & (i0 < 23.0)
        priority[w2] = 2
        logger.info(f'{(keep & w2).sum()} {self.name} stars are marked as priority 2')

        # Priority 3:
        w3 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.8) & (i0 < 22.25)
        priority[w3] = 3
        logger.info(f'{(keep & w3).sum()} {self.name} stars are marked as priority 3')

        # Priority 4:
        w4 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.8) & (i0 < 23.0)
        priority[w4] = 4
        logger.info(f'{(keep & w4).sum()} {self.name} stars are marked as priority 4')

        # Priority 5:
        w5 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.1) & (i0 < 22.25)
        priority[w5] = 5
        logger.info(f'{(keep & w5).sum()} {self.name} stars are marked as priority 5')

        # Priority 6:
        w6 = (priority == -1) & np.isfinite(p_member) & (p_member > 0.1) & (i0 < 23.0)
        priority[w6] = 6
        logger.info(f'{(keep & w6).sum()} {self.name} stars are marked as priority 6')

        # Only keep stars with valid priority
        keep &= (priority >= 0) & (priority <= 9) & (g0 > 19.0) & (i0 > 19.0) & (i0 < 30.0) & (g0 < 30.0)

        catalog.data['priority'] = -1
        catalog.data.loc[keep, 'priority'] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data.loc[keep, 'exp_time'] = exp_time[keep]


    def assign_priorities_old(self, catalog: Catalog, mask=None, isogrid=None):
        """Assign priority classes based on photometry"""

        mask = mask if mask is not None else np.s_[:]

        g0, i0, gi0 = self._get_hsc_dered_mags_colors(catalog, mask)
        clg = catalog.data['clg'][mask]
        cli = catalog.data['cli'][mask]

        priority = np.full(g0.shape, -1, np.int32)
        code = np.full(g0.shape, 0, dtype=np.int32)

        if 'p_member' in catalog.data:
            prob = catalog.data['p_member'][mask]

            top_pri = np.maximum(np.floor((i0 - 21)/(24.5 - 21) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((i0 - 21)/(24.5 - 21) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness

            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)

            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # Bright stars
            w = (i0 <= 21.0) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = (i0 >= 24.0) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 2

            # Possible extended sources, regardless of magnitude
            w = (cli > 0.5) | (clg > 0.5)
            priority[w] = 9
            code[w] = 3

        else:
            raise(NotImplementedError)
        
        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10**((i0-19.0)/2.5)) + 1).astype(int), 10), 1)

        keep = (i0 < 25.0) & (g0 < 24.5) & (priority <= 9) & (code <= 2)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]
        
        catalog.data.loc[(0 <= catalog.data['priority']) & catalog.data['exp_time'].isna(), 'exp_time'] = 1800