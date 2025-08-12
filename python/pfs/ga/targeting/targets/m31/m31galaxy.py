import astropy.units as u

from ...util.args import *
from ...instrument import *
from ...projection import Pointing
from ...data import Catalog, Observation
from ...diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from ...photometry import Photometry, Magnitude, Color
from ...selection import ColorSelection, MagnitudeSelection, LinearSelection
from ...diagram import CMD, CCD, ColorAxis, MagnitudeAxis
from ...photometry import Photometry, Magnitude, Color
from ...config.sample import SampleConfig
from ..galaxy import Galaxy

class M31Galaxy(Galaxy):
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
    
    def get_selection_mask(self, catalog: Catalog, nb=True, blue=False, probcut=None, observed=None, bright=16, faint=23.5):
        raise NotImplementedError()
    
    def get_sample_config(self):
        config = SampleConfig()
        return config
