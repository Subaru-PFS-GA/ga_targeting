import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from pfs.ga.targeting.scripts.pmap.pmapscript import PMapScript
from pfs.ga.targeting.config.pmap import PMapConfig
from pfs.ga.targeting.instrument import *


def load_pmap_config(config_file, format='.json'):
    config = PMapConfig.default()
    config.load(config_file, ignore_collisions=True, format=format)
    return config

def load_observations(field, config):
    r = field.get_text_observation_reader()
    obs = r.read(os.path.expandvars(config.obs_path))
    return obs

def load_simulation(config):
    s = PMapScript()
    s._config = config
    return s._PMapScript__load_simulation()