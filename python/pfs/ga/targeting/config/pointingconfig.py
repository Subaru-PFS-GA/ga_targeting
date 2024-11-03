from datetime import datetime

from ..core import Pointing

from .config import Config

class PointingConfig(Config):
    def __init__(self,
                 ra: float = None,
                 dec: float = None,
                 posang: float = None,
                 obs_time: datetime = None,
                 exp_time: float = None):
        
        self.ra = ra
        self.dec = dec
        self.posang = posang
        self.obs_time = obs_time
        self.exp_time = exp_time

        super().__init__()

    def get_pointing(self):
        return Pointing(self.ra, self.dec, posang=self.posang, obs_time=self.obs_time, exp_time=self.exp_time)
    