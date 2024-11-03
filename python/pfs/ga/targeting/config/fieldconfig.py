from typing import List, Dict
from datetime import datetime

from .config import Config

class FieldConfig(Config):
    def __init__(self,
                 key: str = None,
                 name: str = None,
                 arms: str = None,
                 nvisits: int = None,
                 exp_time: float = None,
                 obs_time: datetime = None):
        
        # Short name for the main object of the field
        self.key = key

        # Full name of the field
        self.name = name

        # List of arms to be used to observe the field. It has no effect on the
        # targeting algorithm, but these values are used to generate the pfsDesign files.
        self.arms = arms

        # Number of visits for each pointing.
        self.nvisits = nvisits

        # Exposure time for each visit, in seconds.
        self.exp_time = exp_time

        # Obsevation time, UTC, for all visits. Although objects move in the focal plane over
        # time, the difference in position is assumed to be negligible for the purposes of targeting.
        # The pfsDesign files will be further tweaked to position the fibers to the correct positions.
        self.obs_time = obs_time

        super().__init__()