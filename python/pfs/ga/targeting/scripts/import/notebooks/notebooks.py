import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

from pfs.ga.targeting.config.netflow import NetflowConfig
from pfs.ga.targeting.scripts.netflow.netflowscript import NetflowScript
from pfs.ga.targeting.projection import WcsProjection, Pointing
from pfs.ga.targeting.diagram import CMD, CCD, FOV, FP, ColorAxis, MagnitudeAxis
from pfs.ga.targeting.diagram import CMD, CCD, FOV, FP, ColorAxis, MagnitudeAxis
from pfs.ga.targeting.photometry import Photometry, Magnitude, Color
from pfs.ga.targeting.instrument import *

def load_netflow_config(config_file, format='.json'):
    config = NetflowConfig.default()
    config.load(config_file, ignore_collisions=True, format=format)
    return config

def print_target_list_config(config):
    for key in config.targets:
        print(key, config.targets[key].prefix, config.targets[key].path)

def load_target_lists(config, output_path, prefix=None):
    target_lists = {}
    for key in config.targets:
        if prefix is None or config.targets[key].prefix in prefix:
            # Try loading from the output directory, then from the source directory
            fn = os.path.expandvars(os.path.join(output_path, f'{config.field.key}_targets_{key}.feather'))
            if not os.path.isfile(fn):
                fn = os.path.expandvars(config.targets[key].path)
            target_lists[key] = NetflowScript.load_target_list(key, config.targets[key], fn)
    return target_lists

def print_photometry(target_lists):
    # List available photometry for each target list
    for k, target_list in target_lists.items():
        print(k)
        for p in target_list.photometry:
            print(' ', p)
            for m in target_list.photometry[p].magnitudes:
                print('    ', m)

def get_unique_photometry(target_lists):
    # Unique photometry across all target lists
    phot = np.unique([p for k in target_lists
                        for p in target_lists[k].photometry.keys()
                        if target_lists[k].photometry is not None])
    
    # Unique magnitudes across all photometry
    mags = {}
    for p in phot:
        mm = []
        for k, target_list in target_lists.items():
            if target_list.photometry is not None and p in target_list.photometry:
                mm.extend([ m for m in list(target_list.photometry[p].magnitudes) ])
        mags[p] = np.unique(mm)

    return phot, mags

def create_fov(center):
    wcs = WcsProjection(center, proj='TAN')
    wfc = SubaruWFC(center)
    fov = FOV(projection=wcs)

    return wcs, wfc, fov

def create_cmd(photometry):
    # Create a color-magnitude diagram using the first two filters of the photometric systems
    cmd = {}
    for p, phot in photometry.items():
        if len(phot.magnitudes) > 1:
            mm = list(phot.magnitudes.keys())
            m1, m2 = phot.magnitudes[mm[0]], phot.magnitudes[mm[1]]
            cmd[p] = CMD([ColorAxis(Color([m1, m2])), MagnitudeAxis(m2)])
    return cmd

def plot_target_list_coordinates(ax, fov, target_list, **kwargs):

    n = len(target_list.data)
    N = 100000
    if n <= N:
        alpha = 1.0
    else:
        alpha = max(0.01, np.exp(-(n - N) / N))

    l = fov.plot_catalog(ax, target_list, size=1, scalex=False, scaley=False, alpha=alpha, **kwargs)

    ax.set_title(target_list.name)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlim(3, -3)
    ax.set_ylim(-3, 3)
    ax.grid(True)

    return l

def plot_pointings(ax, pfi, wcs, wfc, fov, pointings, **kwargs):
    for pointing in pointings:
        # wcs = WcsProjection(pointing, proj='TAN')
        wfc = SubaruWFC(pointing)
        # fov = FOV(projection=wcs)
        pfi.plot_focal_plane(ax, fov, corners=True, projection=wfc, color='red')

def plot_target_list_cmd(ax, target_list, photometry, m1, m2, **kwargs):
    cmd = CMD([ColorAxis(Color([m1, m2])), MagnitudeAxis(m2)])
    
    l = cmd.plot_catalog(ax, target_list, size=1, **kwargs)

    ax.set_title(photometry.name)
    ax.grid(True)

    return l

def plot_target_list(target_list, pfi, center, pointings=None):
    # Color by priority, if available
    if 'priority' in target_list.data:
        color = target_list.data['priority']
    else:
        color = 'black'

    wcs, wfc, fov = create_fov(center)

    cmap = plt.get_cmap('tab10')
    f = plt.figure(figsize=(8, 4), dpi=240)
    gs = GridSpec(1, 2, figure=f)
    axs = [None, None]

    # Add an axis based on gridspec
    axs[0] = f.add_subplot(gs[0], projection=wcs.wcs)
    axs[1] = f.add_subplot(gs[1])

    l = plot_target_list_coordinates(axs[0], fov, target_list,
                                     c=color, cmap=cmap)
    
    if pointings is not None:
        plot_pointings(axs[0], pfi, wcs, wfc, fov, pointings, c='red')

    axs[0].set_xlim(1.5, -1.5)
    axs[0].set_ylim(-1.5, 1.5)

    # Plot a color-magnitude diagram using the first two filters of the first photometric
    # system available in the target list that has at least two filters defined in the data set
    for p, photometry in target_list.photometry.items():
        if len(photometry.magnitudes) > 1:
            mm = list(photometry.magnitudes.keys())
            m1, m2 = mm[0], mm[1]
            plot_target_list_cmd(axs[1], target_list, photometry, photometry.magnitudes[m1], photometry.magnitudes[m2], c=color, cmap=cmap)

            break

    f.colorbar(l, ax=axs[1], orientation='vertical', label='Priority')
    f.suptitle(target_list.name)
    f.show()

def plot_magnitude_dist(target_lists, phot, mags):
    # Plot the magnitude distribution for each target list for each magnitude
    cmap = plt.get_cmap('tab10')

    for p in phot:
        for m in mags[p]:
            f, ax = plt.subplots(1, 1, figsize=(3.5, 2.5), dpi=240)

            for i, (k, target_list) in enumerate(target_lists.items()):
                if target_list.photometry is not None and p in target_list.photometry and m in target_list.photometry[p].magnitudes:
                    mag, mag_err = target_list.get_magnitude(target_list.photometry[p].magnitudes[m])
                    mask = np.isfinite(mag)
                    hist, bins = np.histogram(mag[mask], bins=30)
                    ax.plot(bins[:-1], hist, color=cmap(i), label=target_list.name)
                        
            ax.set_xlabel(f'{p} {m} magnitude')
            ax.set_ylabel('Number of targets')
            ax.legend()

def plot_priority_dist(target_list):
    # Plot the priority distribution for each target list
    f, ax = plt.subplots(1, 1, figsize=(3.5, 2.5), dpi=240)

    hist = np.bincount(target_list.data['priority'])
    ax.bar(np.arange(hist.size), hist)

    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax.set_title(f'Priority class distribution for target list `{target_list.name}`')
    ax.set_xlabel(f'Priority class')
    ax.set_ylabel(f'Target count')

def plot_required_visits_dist(target_list, exp_time):
    f, ax = plt.subplots(1, 1, figsize=(3.5, 2.5), dpi=240)

    hist = np.bincount((np.ceil(target_list.data['exp_time'] / exp_time)).astype(int))
    ax.bar(np.arange(hist.size), hist)

    # ax.set_xlim(0.01, None)

    # Only plot integer values along x-axis
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax.set_title(f'Distribution of required visit (t_exp = {exp_time} s)')
    ax.set_xlabel('Number of required visits')
    ax.set_ylabel('Target count')