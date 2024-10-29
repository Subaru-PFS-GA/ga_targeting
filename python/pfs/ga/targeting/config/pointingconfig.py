from datetime import datetime

from ..core import Pointing

from .config import Config

class PointingConfig(Config):
    def __init__(self):
        
        self.ra = None
        self.dec = None
        self.posang = None
        self.obs_time = None
        self.exp_time = None

        super().__init__()

    def get_pointing(self):
        return Pointing(self.ra, self.dec, posang=self.posang, obs_time=self.obs_time, exp_time=self.exp_time)
    