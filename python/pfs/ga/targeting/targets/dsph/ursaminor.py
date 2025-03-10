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
from ...config.pmap import PMapConfig
from ...config.sample import SampleConfig
from .dsphgalaxy import DSphGalaxy

prioritize_DEIMOS = True

class UrsaMinor(DSphGalaxy):
    def __init__(self):
        ID = 'umi'
        name = 'Ursa Minor'
        # pos = [ '15h 09m 08.5s', '+67d 13m 21s' ]                   # Simbad, Carlsten et al. (2021)
        # pos = [ 227.29725, 67.21436111 ] * u.deg                    # Evan
        
        # Radial profile fit parameters
        pos = [ 227.27229069,  67.23449905 ] * u.deg                  # Profile fit, Dobos
        white = (np.array([-0.0433267, -0.02317608]), np.array([[[-4.37100736, -3.55293437], [-9.83032195, 12.09378084]]]))
        king = (45.302893030209475, 1084.807251671976, 1.3825051586583355, 7.252224797182235)

        rad = 120 * u.arcmin                                          # Approximate field radius
        
        DM, DM_err = 18.9, 0.2
        pm = [ -0.119, 0.072 ] * u.mas / u.yr                         # Pace et al. (2022)
        pm_err = [ 0.005, 0.005 ] * u.mas / u.yr
        RV, RV_err = (-274.0, 1.0) * u.kilometer / u.second

        # Naive pointings
        # ra0 = [ 229.2, 226.0, 225.2, 228.5, 228.2, 226.3, 226.0, 228.0 ] * u.deg
        # dec0 = [ 67.90, 67.80, 66.55, 66.60, 67.5, 67.5, 66.9, 66.955 ] * u.deg
        # pa = [ 0, 0, 0, 0, 0, 0, 0, 0] * u.deg
        # pointings = {
        #     SubaruPFI: [ Pointing((ra, dec), posang=pa) for ra, dec, pa in zip(ra0, dec0, pa) ]
        # }

        # Generate the poitings algorithmically
        pointings = {
            SubaruPFI: [
                Pointing.from_relative_pos(pos, sep=0.35, dir=50, posang=-20),
                Pointing.from_relative_pos(pos, sep=-0.35, dir=50, posang=-20),

                # 6-fold overlap in the center
                # Pointing.from_relative_pos(pos, sep=0.6, dir=50, posang=-20),
                # Pointing.from_relative_pos(pos, sep=-0.6, dir=50, posang=-20),

                # Rotated outer pointings along the minor axis
                # Pointing.from_relative_pos(pos, sep=0.6, dir=50, posang=-50),
                # Pointing.from_relative_pos(pos, sep=-0.6, dir=50, posang=-50),

                # Rotated outer pointings along the minor axis
                Pointing.from_relative_pos(pos, sep=0.85, dir=50, posang=-50),
                Pointing.from_relative_pos(pos, sep=-0.85, dir=50, posang=-50),

                Pointing.from_relative_pos(pos, sep=0.45, dir=140, posang=-20),
                Pointing.from_relative_pos(pos, sep=-0.45, dir=140, posang=-20),

                Pointing.from_relative_pos(pos, sep=1.05, dir=140, posang=-20),
                Pointing.from_relative_pos(pos, sep=-1.05, dir=140, posang=-20),
            ]
        }

        super().__init__(ID, name,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings)
        
        # CMD and CCD definitions with color and magnitude limits for Ursa Minor

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
            center = PointingConfig.from_pointing(self.get_center()),
            arms = 'bmn',
            nvisits = 1,
            exp_time = 30 * 60.,        # 3 hr total
            obs_time = datetime(2025, 5, 25, 0, 0, 0) + timedelta(hours=10),
            resolution = 'm',
        )

        config.pointings = [ PointingConfig(p.ra, p.dec, p.posang) for p in self.get_pointings(SubaruPFI) ]

        return config

    def get_pmap_config(self):
        config = PMapConfig(
            cut_nb = True,
            keep_blue = True,
            extents = [[0.1, 2.0], [17.0, 23.5]],
            merge_list = [np.s_[:10], np.s_[10:]]
        )

        return config

    def get_sample_config(self):
        config = SampleConfig(

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
            # 'i_hsc': 'i2_hsc',
        }

        return filter_map
    
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
                & ColorSelection(ccd.axes[0], None, 1.65).apply(catalog, observed=observed)
                
                | LinearSelection(ccd.axes, [-0.25, 1.0], -0.15, None).apply(catalog, observed=observed)
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

            top_pri = np.maximum(np.floor((g0 - 17)/(23.0 - 17) * 4).astype(int), 0) # top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((g0 - 17)/(23.0 - 17) * 4).astype(int) + 4, 4) # bot pri goes from 4-8 based on brightness
          
            w = ~np.isnan(prob)
            priority[w] = np.minimum(np.maximum(bot_pri[w] - np.rint(prob[w] * (bot_pri[w] - top_pri[w])).astype(int), 0), 8)
            
            # Everything without membership probability
            w = np.isnan(prob) | (prob == 0.0)
            priority[w] = 9
            code[w] = 0
            
            # Blue Horizontal Branch
            w = (g0 > 19.4) & (g0 < 20.2) & (gi0 > -0.5) & (gi0 < 0.2) & (cli <= 0.5) & (clg < 0.5)
            priority[w] = 6
            code[w] = 0
            
            # Blue Stragglers
            x1,x2,x3,x4 = -0.7, -0.2, -0.2, 0.3
            y1,y2,y3,y4 = 20.8, 20.8, 23.0, 23.0
            w = (priority == 9) & (g0 > y1) & (g0 < y3) & ((g0-y1) < ((y3-y1)/(x3-x1)*(gi0-x1))) & ((g0-y2) > ((y4-y2)/(x4-x2)*(gi0-x2)))
            priority[w] = 7
            code[w] = 0
            
            # Very bright stars, this does nothing because there aren't any of those
            w = (i0 <= 16) & (cli <= 0.5) & (clg <= 0.5)
            priority[w] = 9
            code[w] = 1
            
            # Very faint stars with lowest priority
            w = ((g0 >= 23.0) | (i0 >= 23.0)) & (cli <= 0.5) & (clg <= 0.5)
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

        
        #If a target has been previously observed by DEIMOS (Kirby et al. 2010), then set its priority to 0.
        if prioritize_DEIMOS:
            import os
            from astropy.coordinates import SkyCoord
            from astropy.io import fits

            fn = os.path.expandvars('$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/umi_moogify_member.fits.gz')
            hdul = fits.open(fn)
            deimos = hdul[1].data
            c_deimos = SkyCoord(deimos['RA']*u.degree, deimos['DEC']*u.degree)
            c_hsc = catalog.get_skycoords() #SkyCoord(catalog.data['RA']*u.degree, catalog.data['Dec']*u.degree)
            idx, d2d, _ = c_deimos.match_to_catalog_sky(c_hsc)
            priority[idx] = 0
            code[idx] = 0


        exp_time = 1800 * np.maximum(np.minimum(np.rint(5 * np.sqrt(10**((i0-19.0)/2.5)) + 1).astype(int), 6), 1)

        keep = (g0 < 23) & (priority <= 9) & (code == 0)

        catalog.data['priority'] = -1
        catalog.data['priority'][mask][keep] = priority[keep]

        catalog.data['exp_time'] = np.nan
        catalog.data['exp_time'][mask][keep] = exp_time[keep]