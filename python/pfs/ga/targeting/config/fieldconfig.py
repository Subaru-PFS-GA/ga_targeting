from typing import List, Dict
from datetime import datetime

from .config import Config
from .pointingconfig import PointingConfig

class FieldConfig(Config):
    def __init__(self,
                 pointings: List[PointingConfig] = None):
        
        self.key = None
        self.name = None
        self.definition = None
        self.pointings = pointings
        self.arms = None
        self.nvisits = None
        self.exp_time = None
        self.obs_time = None

        super().__init__()