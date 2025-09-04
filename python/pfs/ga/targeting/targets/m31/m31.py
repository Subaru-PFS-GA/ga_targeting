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
        # group_name, field_name,  ra,             dec, posang, stage, priority
        ('m31_E0', 'diskE_1',    '00:51:50.50',  '+42:39:59.8', 170.0, 0, 1),
        ('m31_E0', 'diskE_2',    '00:50:39.13',  '+41:32:16.9', 170.0, 0, 2),
        ('m31_E0', 'haloE_1',    '00:56:29.68',  '+41:54:14.8', 170.0, 0, 5),

        ('m31_W0', 'diskW_2',    '00:40:30.97',  '+42:45:33.3', 180.0, 0, 6),
        ('m31_W0', 'diskW_3',    '00:34:52.69',  '+42:15:58.8', 180.0, 0, 7),
        ('m31_W0', 'haloW_4',    '00:35:20.12',  '+43:24:47.6', 180.0, 0, 8),
        
        ('m31_GSS0', 'GSS_1',      '00:44:27.95',  '+39:41:48.1', 170.0, 0, 4),

        ('m31_NWS0', 'NWstream_5', '00:10:34.58',  '+46:50:16.5', 180.0, 0, 3),

        ##
        # TODO: make sure the position angles are updated below
        
        ('none', 'GSS_2',      '00:38:59.79',  '+39:14:15.7', 110.0, 0, 0),
        ('none', 'GSS_3',      '00:48:23.81',  '+39:16:42.8', 110.0, 0, 0),
        ('none', 'GSS_4',      '00:42:49.51',  '+38:53:14.9', 110.0, 0, 0),
        ('none', 'GSS_5',      '00:52:39.04',  '+38:30:07.3', 110.0, 0, 0),
        ('none', 'GSS_6',      '00:47:33.81',  '+38:09:19.1', 110.0, 0, 0),

        ('none', 'diskW_1',    '00:41:04.96',  '+43:54:16.6', 120.0, 0, 0),
        ('none', 'diskW_4',    '00:34:26.15',  '+41:07:09.7', 120.0, 0, 0),

        ('none', 'haloW_1',    '00:37:29.04',  '+46:23:47.7', 120.0, 0, 0),
        ('none', 'haloW_2',    '00:36:56.68',  '+45:14:48.3', 120.0, 0, 0),
        ('none', 'haloW_3',    '00:36:11.59',  '+44:12:09.2', 120.0, 0, 0),
        ('none', 'haloW_5',    '00:30:15.25',  '+44:04:51.9', 120.0, 0, 0),
        ('none', 'haloW_6',    '00:29:40.85',  '+42:54:16.5', 120.0, 0, 0),
        ('none', 'haloW_7',    '00:29:19.71',  '+41:45:24.0', 120.0, 0, 0),
        ('none', 'haloW_8',    '00:24:38.32',  '+44:40:35.0', 120.0, 0, 0),
        ('none', 'haloW_9',    '00:24:22.57',  '+43:31:40.3', 120.0, 0, 0),
        ('none', 'haloW_10',   '00:24:38.91',  '+42:26:52.2', 120.0, 0, 0),
        ('none', 'haloW_11',   '00:18:57.73',  '+44:08:07.7', 120.0, 0, 0),
        ('none', 'haloW_12',   '00:18:48.36',  '+42:59:12.3', 120.0, 0, 0),
        
        ('none', 'diskE_3',    '00:50:01.64',  '+40:24:31.2', 110.0, 0, 0),

        ('none', 'haloE_2',    '00:55:14.42',  '+40:46:43.3', 110.0, 0, 0),
        ('none', 'haloE_3',    '00:54:01.84',  '+39:39:08.5', 110.0, 0, 0),
        ('none', 'haloE_4',    '01:01:01.75',  '+41:04:59.3', 110.0, 0, 0),
        ('none', 'haloE_5',    '00:59:49.28',  '+40:03:07.7', 110.0, 0, 0),
        ('none', 'haloE_6',    '00:58:27.46',  '+38:53:07.1', 110.0, 0, 0),
        
        ('none', 'NWstream_1', '00:13:23.24',  '+43:34:41.9', 120.0, 0, 0),
        ('none', 'NWstream_2', '00:08:46.61',  '+44:12:07.3', 120.0, 0, 0),
        ('none', 'NWstream_3', '00:13:26.22',  '+44:43:36.5', 120.0, 0, 0),
        ('none', 'NWstream_4', '00:10:09.21',  '+45:41:58.7', 120.0, 0, 0),
        
    ]

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
                exp_time=30*60, nvisits=10,
                label=fl)
        
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
        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10 ** ((i0 - 19.0) / 2.5)) + 1).astype(int), 10), 1)

        # Priorities
        priority = np.full_like(p_member, -1, np.int32)

        # Everything without membership probability
        w9 = (np.isnan(p_member) | (p_member <= 0.1)) & (i0 < 23.0) & (i0 > 20.0)
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