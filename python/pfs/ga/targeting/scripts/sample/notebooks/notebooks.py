import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from pfs.ga.targeting.scripts.sample.samplescript import SampleScript
from pfs.ga.targeting.config.sample import SampleConfig
from pfs.ga.targeting.instrument import *

def load_sample_config(config_file, format='.json'):
    config = SampleConfig.default()
    config.load(config_file, ignore_collisions=True, format=format)
    return config

def load_observations(field, config):
    r = field.get_text_observation_reader()
    obs = r.read(os.path.expandvars(config.obs_path))
    return obs

def plot_sample(field, sample, cmd, ccd, pfi, fov, wcs,
                mask=None,
                color_by=None, cmap='tab10', label='target priority', title='HSC targets',
                **kwargs):
    
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    f = plt.figure(figsize=(8, 3), dpi=240)
    gs = f.add_gridspec(1, 3, width_ratios=[2, 2, 3], wspace=0.4)

    axs = [f.add_subplot(gs[0, 0]), f.add_subplot(gs[0, 1]), f.add_subplot(gs[0, 2], projection=wcs.wcs)]

    if color_by is not None:
        c = sample.data[color_by][mask if mask is not None else slice(None)]
    else:
        c = None

    cmd.plot_observation(axs[0], sample, c=c, observed=True, mask=mask, cmap=cmap, **kwargs)
    ccd.plot_observation(axs[1], sample, c=c, observed=True, mask=mask, cmap=cmap, **kwargs)
    l = fov.plot_observation(axs[2], sample, c=c, observed=True, mask=mask, cmap=cmap, **kwargs)
            
    pp = field.get_pointings(SubaruPFI)
    for p in pp:
        pfi.plot_focal_plane(axs[2], fov, corners=True, projection=SubaruWFC(p))

    axs[2].set_xlim(0.75, -2.75)
    axs[2].set_ylim(-2.75, 0.75)

    # Put a colorbar under all three subplots
    if c is not None:
        f.colorbar(l, ax=axs[2], label=label)

    f.suptitle(title)
