import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from pfs.ga.targeting.scripts.priority.priorityscript import PriorityScript
from pfs.ga.targeting.config.priority import PriorityConfig
from pfs.ga.targeting.instrument import *

def load_priority_config(config_file, format='.json'):
    config = PriorityConfig.default()
    config.load(config_file, ignore_collisions=True, format=format)
    return config

def load_observations(field, config):
    r = field.get_text_observation_reader()
    obs = r.read(os.path.expandvars(config.obs_path))
    return obs

