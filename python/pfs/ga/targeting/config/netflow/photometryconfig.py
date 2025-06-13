from typing import Dict

from ..config import Config

class PhotometryConfig(Config):
    def __init__(self):

        # Defines the magnitudes. See example.py for more details
        self.filters = None

        # Defines the filter bands. See example.py for more details
        self.bands = None

        # Defines magnitude limits. See example.py for more details
        self.limits = None

        # Mapping between filter names in the input data files and the
        # filter names in the exported target lists.
        self.filter_map = None

        super().__init__()
