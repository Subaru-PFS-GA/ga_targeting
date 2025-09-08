from typing import List, Dict
from datetime import datetime

from pfs.ga.common.config import Config

from .pointingconfig import PointingConfig

class FieldConfig(Config):
    def __init__(self,
                 key: str = None,
                 sector: str = None,
                 name: str = None,
                 id_prefix: int = None,
                 center: PointingConfig = None,
                 arms: str = None,
                 nvisits: int = None,
                 nrepeats: int = None,
                 exp_time: float = None,
                 obs_time: datetime = None,
                 resolution: str = None):
        
        # Short name for the main object of the field
        self.key = key

        # Sector within the field, used for M31
        self.sector = sector

        # Full name of the field
        self.name = name

        # A bit prefix to combine with target_idx to generate
        # unique database identifiers
        self.id_prefix = id_prefix

        # Field center
        self.center = center

        # List of arms to be used to observe the field. It has no effect on the
        # targeting algorithm, but these values are used to generate the pfsDesign files.
        self.arms = arms

        # Number of visits for each pointing.
        self.nvisits = nvisits

        # Number of repeats for each visit with the same design
        # Assumed 1 of set to None
        self.nrepeats = nrepeats

        # Exposure time for each visit, in seconds.
        self.exp_time = exp_time

        # Obsevation time, UTC, for all visits. Although objects move in the focal plane over
        # time, the difference in position is assumed to be negligible for the purposes of targeting.
        # The pfsDesign files will be further tweaked to position the fibers to the correct positions.
        self.obs_time = obs_time

        # Spectral resolution, "low" or "medium" in the red arm
        self.resolution = resolution

        super().__init__()