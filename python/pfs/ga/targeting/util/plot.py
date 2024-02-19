import numpy as np

def find_plot_limits(data, vmin=None, vmax=None):
    mask = np.isfinite(data)
    if mask.sum() > 0:
        vmin = vmin if vmin is not None else np.min(data[mask])
        vmax = vmax if vmax is not None else np.max(data[mask])

    return vmin, vmax

def get_plot_normalized_color(cmap, data, vmin=None, vmax=None):
    if vmin is None or vmax is None:
        return cmap(data)
    else:
        return cmap((data - vmin) / (vmax - vmin))