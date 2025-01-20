from typing import Dict
from .config import Config

class PhotometryConfig(Config):
    def __init__(self):

        # Defines the magnitudes. See example.py for more details
        self.filters = None

        # Defines the filter bands. See example.py for more details
        self.bands = None

        # Defines magnitude limits. See example.py for more details
        self.limits = None

        super().__init__()