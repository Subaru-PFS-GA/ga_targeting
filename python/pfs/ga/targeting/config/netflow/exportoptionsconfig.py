from typing import Dict
import numpy as np

from ..config import Config, Lambda

class ExportOptionsConfig(Config):

    def __init__(self):
        pass

        # Mapping between filters and bands
        self.filter_bands = None

    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        config = ExportOptionsConfig()
        config.filter_bands = {
            'g_hsc': 'g',
            'g_ps1': 'g',
            'g_sdss': 'g',
            'bp_gaia': 'g',

            'r_hsc': 'r',
            'r2_hsc': 'r',
            'r_old_hsc': 'r',
            'r_ps1': 'r',
            'r_sdss': 'r',
            'g_gaia': 'r',

            'i_hsc': 'i',
            'i2_hsc': 'i',
            'i_old_hsc': 'i',
            'i_ps1': 'i',
            'i_sdss': 'i',
            'rp_gaia': 'i',
                        
            'z_hsc': 'z',
            'z_ps1': 'z',
            'z_sdss': 'z',

            'y_hsc': 'y',
            'y_ps1': 'y',
        }
        return config