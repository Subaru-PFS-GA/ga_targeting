from datetime import datetime

from pfs.ga.common.projection import Pointing
from pfs.ga.common.config import Config

class PointingConfig(Config):
    def __init__(self,
                 ra: float = None,
                 dec: float = None,
                 posang: float = None,
                 obs_time: datetime = None,
                 exp_time: float = None,
                 priority: int = None,
                 stage: int = None):
        
        self.ra = ra
        self.dec = dec
        self.posang = posang
        self.obs_time = obs_time
        self.exp_time = exp_time
        self.priority = priority
        self.stage = stage

        super().__init__()

    @staticmethod
    def from_pointing(pointing: Pointing):
        return PointingConfig(
            ra=pointing.ra,
            dec=pointing.dec,
            posang=pointing.posang,
            obs_time = pointing.obs_time.value if pointing.obs_time is not None else None,
            exp_time = pointing.exp_time.to_value('s') if pointing.exp_time is not None else None,
            priority = pointing.priority,
            stage = pointing.stage
        )

    def get_pointing(self, obs_time=None, exp_time=None):
        obs_time = obs_time if obs_time is not None else self.obs_time
        exp_time = exp_time if exp_time is not None else self.exp_time

        return Pointing(self.ra, self.dec, posang=self.posang, obs_time=obs_time, exp_time=exp_time,
                        priority=self.priority, stage=self.stage)
    