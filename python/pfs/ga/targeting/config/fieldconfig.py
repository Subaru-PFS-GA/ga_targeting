from typing import List, Dict
from datetime import datetime

from .config import Config

class FieldConfig(Config):
    def __init__(self):
        
        self.key = None
        self.name = None
        self.arms = None
        self.nvisits = None
        self.exp_time = None
        self.obs_time = None

        super().__init__()