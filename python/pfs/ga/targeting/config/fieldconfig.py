from typing import List, Dict
from datetime import datetime

from .config import Config

class FieldConfig(Config):
    def __init__(self):
        
        # Short name for the main object of the field
        self.key = None

        # Full name of the field
        self.name = None

        # List of arms to be used to observe the field. It has no effect on the
        # targeting algorithm, but these values are used to generate the pfsDesign files.
        self.arms = None

        # Number of visits for each pointing.
        self.nvisits = None

        # Exposure time for each visit, in seconds.
        self.exp_time = None

        # Obsevation time, UTC, for all visits. Although objects move in the focal plane over
        # time, the difference in position is assumed to be negligible for the purposes of targeting.
        # The pfsDesign files will be further tweaked to position the fibers to the correct positions.
        self.obs_time = None

        super().__init__()