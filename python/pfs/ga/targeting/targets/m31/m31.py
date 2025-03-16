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
from .m31galaxy import M31Galaxy

class M31(M31Galaxy):
    def __init__(self):
        ID = 'm31'
        name = 'M31'

        pos = [ '00h 42m 44.33s', '+41d 16m 07.5s' ]
        rad = 6 * u.deg
        DM, DM_err = 24.407, 0.032   #Siyang Li et al. (2021, ApJ, 920, 84)
        pm = [ -0.0533, -0.0104 ] * u.mas / u.yr    #Sohn et al. (2012)
        pm_err = [ 0.0246, 0.0244 ] * u.mas / u.yr
        RV, RV_err = (-287.2, 8.0) * u.kilometer / u.second

        rastr = [ '00:34:54',
                '00:28:56',
                '00:22:49',
                '00:41:43',
                '00:35:47',
                '00:29:42',
                '00:23:26',
                '00:42:45',
                '00:36:42',
                '00:37:40',
                '00:34:04',
                '00:28:13',
                '00:22:13',
                '00:50:18',
                '00:55:44',
                '01:01:02',
                '00:51:28',
                '00:56:59',
                '01:02:22',
                '00:52:41',
                '00:58:18',
                '01:03:46',
                '00:43:43',
                '00:49:12',
                '00:54:32',
                '00:37:13',
                '00:42:44',
                '00:48:06',
                '00:16:31.7',
                '00:10:04.7',
                '00:10:24.5',
                '00:16:04.2',
                '00:09:45.8',
                '00:06:10',
                '00:52:10',
                '00:48:10',
                '00:47:06',
                '00:47:05',
                '00:39:30',
                '00:37:40',
                '00:32:40',
                '00:31:10' ]
        decstr = [ '+42:22:57',
                '+43:11:01',
                '+43:57:54',
                '+42:54:25',
                '+43:43:47',
                '+44:32:01',
                '+45:19:01',
                '+44:15:02',
                '+45:04:35',
                '+46:25:20',
                '+41:02:05',
                '+41:50:00',
                '+42:36:45',
                '+40:07:25',
                '+39:15:27',
                '+38:22:34',
                '+41:27:45',
                '+40:35:34',
                '+39:42:28',
                '+42:48:01',
                '+41:55:37',
                '+41:02:17',
                '+39:37:51',
                '+38:47:03',
                '+37:55:17',
                '+39:06:52',
                '+38:17:16',
                '+37:26:39',
                '+44:43:30',
                '+45:27:47',
                '+46:49:07',
                '+43:22:15',
                '+44:06:28',
                '+47:35:07',
                '+44:25:07',
                '+43:30:07',
                '+42:16:07',
                '+40:40:07',
                '+41:50:07',
                '+40:05:04',
                '+40:05:04',
                '+39:05:04' ]
        pa0 = [0]*len(rastr)

        #pointings = {
        #    SubaruPFI: [ Pointing((Angle(ra, unit=u.deg)*15, Angle(dec, unit=u.deg)), posang=pa*u.deg) for ra, dec, pa in zip(rastr, decstr, pa0) ]
        #}

        i = 10
        # just do one pointing for now
        rastr = rastr[i]
        decstr = decstr[i]
        pa = pa0[i]
        pointings = {
            SubaruPFI: [ Pointing((Angle(rastr, unit=u.deg)*15, Angle(decstr, unit=u.deg)), posang=pa*u.deg) ]
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

    def get_netflow_config(self):
        config = NetflowConfig.default()

        config.field = FieldConfig(
            key = self.ID,
            name = self.name,
            arms = 'bmn',
            nvisits = 1,
            exp_time = 30 * 60.,        # 3 hr total
            obs_time = datetime(2025, 1, 24, 20, 0, 0) + timedelta(hours=10),
        )

        config.pointings = [ PointingConfig(p.ra, p.dec, p.posang) for p in self.get_pointings(SubaruPFI) ]

        return config
    
    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=21.5, faint=23.5):
        """Return true for objects within sharp magnitude cuts."""

        # TODO: add Keyi's cut
        
        cmd = self.__hsc_cmd
        ccd = self.__hsc_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], 0.95, 2.5).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            mask &= (
                ColorSelection(ccd.axes[0], 0.12, 0.5).apply(catalog, observed=observed)

                | ColorSelection(ccd.axes[1], -0.05, None).apply(catalog, observed=observed)
                & ColorSelection(ccd.axes[0], None, 2.5).apply(catalog, observed=observed)
                
                #| LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed)
                & LinearSelection(ccd.axes, [1.0, 7.0], 1.8, None).apply(catalog, observed=observed)
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

            top_pri = np.maximum(np.floor((i0 - 20)/(23.0 - 20) * 8).astype(int) - 7, -7) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((i0 - 20)/(23.0 - 20) * 6).astype(int) + 3, 3) # bot pri goes from 3-8 based on brightness
          
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w]**(1/10) * (bot_pri[w] - top_pri[w])).astype(int), 0), 9)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # Very bright stars, this does nothing because there aren't any of those
            w = (i0 <= 16) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = (i0 >= 23.0) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 2

            # Possible extended sources, regardless of magnitude
            w = (cli > 0.5) | (clg > 0.5)
            priority[w] = 9
            code[w] = 3

        else:
            raise(NotImplementedError)
        
        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10**((i0-19.0)/2.5)) + 1).astype(int), 10), 1)

        keep = (g0 < 23) & (priority <= 9) & (code == 0)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]
        
        catalog.data.loc[(0 <= catalog.data['priority']) & catalog.data['exp_time'].isna(), 'exp_time'] = 1800