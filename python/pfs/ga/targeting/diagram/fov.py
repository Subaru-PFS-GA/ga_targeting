from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes.core import WCSAxes

from . import styles
from .spatialdiagram import SpatialDiagram
from .radecaxis import RaDecAxis

class FOV(SpatialDiagram):
    """
    A plot of the field of view of some instrument assuming a valid WCS projection
    defined in the form of WCSAxes or WCSAxesSubplot.
    """

    def __init__(self, projection=None, orig=None):
        ra = RaDecAxis('RA', invert=True)
        dec = RaDecAxis('Dec')
        axes = [ ra, dec ]
        super().__init__(axes, projection=projection, orig=orig)

        self._validate()

    def _validate(self):
        pass

    def _get_native_frame(self, native_frame=None):
        return native_frame if native_frame is not None else 'world'
    
    def apply(self, ax: plt.Axes):
        super().apply(ax)
        ax.set_aspect('equal', adjustable='datalim')

    def plot_radial_profile(self, ax: plt.Axes, profile, R=1, **kwargs):
        # TODO: add default style
        ell = profile.get_ellipse(R)
        return self.plot(ax, ell, native_frame='world', **kwargs)

    def plot_instrument(self, ax: plt.Axes, instrument, **kwargs):
        instrument.plot_field_of_view(ax, self, **kwargs)