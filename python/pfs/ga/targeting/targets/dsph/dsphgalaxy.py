from datetime import datetime, timedelta
import astropy.units as u

from pfs.ga.common.util.args import *
from pfs.ga.common.diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from pfs.ga.common.photometry import Photometry, Magnitude, Color

from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...config.sample import SampleConfig
from ...config.netflow import NetflowConfig, FieldConfig, PointingConfig
from ..galaxy import Galaxy

class DSphGalaxy(Galaxy):
    def __init__(self,
                 ID, name, id_prefix,
                 pos, rad=None,
                 DM=None, DM_err=None,
                 pm=None, pm_err=None,
                 RV=None, RV_err=None,
                 pointings=None,
                 **kwargs):
        
        super().__init__(ID, name, id_prefix,
                         pos, rad=rad,
                         DM=DM, DM_err=DM_err,
                         pm=pm, pm_err=pm_err,
                         RV=RV, RV_err=RV_err,
                         pointings=pointings,
                         **kwargs)

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

    def get_text_observation_reader(self, instrument=SubaruHSC):
        if instrument == SubaruHSC:
            return SubaruHSC.text_observation_reader(
                mags=['i', 'g', 'nb515'], ext=['g', 'i', 'nb515'])
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
            cmd = self._hsc_cmd
        elif instrument == Gaia:
            cmd = self._gaia_cmd
        else:
            raise NotImplementedError()
            
        return cmd
    
    def get_ccd(self, instrument=SubaruHSC):
        """
        Return the definition of the color-color diagram for the given instrument.
        """
        if instrument == SubaruHSC:
            ccd = self._hsc_ccd
        else:
            raise NotImplementedError()

        return ccd
    
    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        raise NotImplementedError()

    def get_pmap_config(self):
        raise NotImplementedError()

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