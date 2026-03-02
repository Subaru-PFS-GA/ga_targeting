from datetime import datetime, timedelta
import astropy.units as u

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...io import ObservationSerializer
from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...config.sample import SampleConfig
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ...config.pmap import PMapConfig
from ..field import Field

class GC(Field):
    def __init__(self,
                 ID, name, id_prefix,
                 pos, rad=None,
                 DM=None, DM_err=None,
                 pm=None, pm_err=None,
                 RV=None, RV_err=None,
                 pointings=None,
                 **kwargs):
        
        super().__init__(ID, name, id_prefix, pos, **kwargs)

        # Bounding radius
        self.__rad = normalize_angle(rad, u.arcmin)

        # Telescope pointings with position and position angle
        self.__pointings = pointings

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

        ps1 = PS1.photometry()
        self._ps1_cmd = CMD([
            ColorAxis(Color([ps1.magnitudes['g'], ps1.magnitudes['r']]), limits=(-0.5, 1.5)),
            MagnitudeAxis(ps1.magnitudes['g'], limits=(14.0, 22.0))
        ])

    def get_text_observation_reader(self, instrument=PS1):
        if instrument == PS1:
            # This is the reader config for the HSC × PS1 × GAIA file from Kota

            reader = ObservationSerializer(format='.csv')
            reader.append_photometry(PS1.photometry())

            reader.column_names = [
                'RA', 'Dec', 'objid',
                'obs_ps1_g', 'err_ps1_g',
                'obs_ps1_r','err_ps1_r',
                'rKmag',
                'obs_ps1_i','err_ps1_i',
                'obs_ps1_z','err_ps1_z',
                'obs_ps1_y', 'err_ps1_y',
                'gPSFf',
                'rPSFf',
                'iPSFf',
                'zPSFf',
                'Epoch',
                'coord',
                'separation',
                'ext_ps1_g',
                'ext_ps1_r',
                'ext_ps1_i',
                'ext_ps1_z',
                'ext_ps1_y'
            ]

            reader.kwargs = dict(
                delimiter=',',
                skiprows = 1
            )

            return reader
        else:
            raise NotImplementedError()

    def get_photometry(self, instrument=PS1):
        """
        Return the photometric system for the given instrument.
        """
        return instrument.photometry()

    def get_cmd(self, index=0, instrument=PS1):
        """
        Return the definition of the color-magnitude diagram for the given instrument.
        """

        if instrument == PS1:
            if index is not None and isinstance(self._ps1_cmd, Iterable):
                cmd = self._ps1_cmd[index]
            else:
                cmd = self._ps1_cmd
        else:
            raise NotImplementedError()
            
        return cmd

    def get_ccd(self, instrument=PS1):
        return None

    def get_pmap_config(self):
        config = PMapConfig(
            cut_nb = False,
            keep_blue = False,
            extents = [[-0.1, 1.0], [14.0, 22.0]],
            merge_list = [np.s_[:10], np.s_[10:]]
        )

        return config

    def get_spatial_selection_mask(self, catalog: Catalog, radius=None):
        return True
    
    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        raise NotImplementedError()
    
    def get_sample_config(self):
        config = SampleConfig()
        return config
    
    def get_selection_mask(self, catalog: Catalog, nb=False, blue=False, probcut=None, observed=None, bright=14.0, faint=22.0):
        """Return true for objects within sharp magnitude cuts."""

        cmd = self._ps1_cmd
        # ccd = self._ps1_ccd

        # Broadband colors
        mask = ColorSelection(cmd.axes[0], -0.1, 1.0).apply(catalog, observed=observed)

        # Narrow band
        if nb:
            raise NotImplementedError()

        # Probability-based cut (map) - nonzero membership probability
        if probcut is not None:
            mask &= probcut.apply(catalog, observed=observed)

        # Allow blue
        if blue:
            raise NotImplementedError

        # Always impose faint and bright magnitude cuts
        mask &= MagnitudeSelection(cmd.axes[1], bright, faint).apply(catalog, observed=observed)

        return mask

    def get_netflow_config(self):
        config = NetflowConfig.default()

        config.field = self.get_field_config()
        config.pointings = [ PointingConfig.from_pointing(p) for p in self.get_pointings(SubaruPFI) ]

        return config

    def get_field_config(self):
        """
        Return the default field configuration for dSph galaxies.
        """

        # TODO: it doesn't check visibility at abs_time

        return FieldConfig(
            key = self.ID,
            name = self.name,
            id_prefix = self.id_prefix,
            center = PointingConfig.from_pointing(self.get_center()),
            arms = 'bmn',
            nvisits = 1,
            exp_time = 30 * 60.,        # 3 hr total
            obs_time = datetime(2025, 5, 25, 0, 0, 0) + timedelta(hours=10),
            resolution = 'm',
        )