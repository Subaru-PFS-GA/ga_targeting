import pandas as pd

from ..core import Pointing
from ..instrument import SubaruPFI
from ..targets.dsph import GALAXIES as DSPH_FIELDS
from ..targets.m31 import M31_FIELDS
from .script import Script

from ..setup_logger import logger

class TargetingScript(Script):
    """
    Base class that implements common function for targeting scripts
    """

    def __init__(self):
        super().__init__()

        self._nvisits = None
        self._field = None
        self._config = None

    def __get_config(self):
        return self._config
    
    config = property(__get_config)

    def _add_args(self):
        super()._add_args()

    def _create_instrument(self):
        return SubaruPFI(instrument_options=self._config.instrument_options)
    
    def _generate_pointings(self):
        # The list of pointings can be defined in the config file, which then
        # override the pointings defined for the field in the library source code.

        if self._config.pointings is not None:
            pp = [ p.get_pointing() for p in self._config.pointings ]
        elif self._field is not None:
            # Load from the class
            pp = self._field.get_pointings(SubaruPFI)
        else:
            raise NotImplementedError()
        
        # The number of visits can be overridden from the command-line argument,
        # see __init_from_args()
        nvisits = self._config.field.nvisits
        exp_time = self._config.field.exp_time
        obs_time = self._config.field.obs_time
        
        # Create new pointing objects with obs_time, exp_time etc.
        pointings = []
        for p in pp:   
            pointings.append(Pointing(
                p.ra, p.dec,
                posang = p.posang,
                obs_time = obs_time,
                exp_time = exp_time,
                nvisits = nvisits))
            
        return pointings

    def _save_assignments_all(self, assignments_all):
        # This dataframe contains columns with dtype `object`` that contain the list of
        # magnitudes, fluxes, filter names, etc.
        
        path = self._get_assignments_all_path()
        logger.info(f'Saving assignments joined with the input catalogs to `{path}`.')
        assignments_all.to_feather(path)

    def _load_assignments_all(self):
        path = self._get_assignments_all_path()
        logger.info(f'Loading assignments joined with the input catalogs from `{path}`.')
        return pd.read_feather(path)