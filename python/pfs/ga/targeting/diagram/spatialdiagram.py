from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes.core import WCSAxes

from ..util import normalize_coords, denormalize_coords
from ..data import Observation, Simulation
from .diagram import Diagram

class SpatialDiagram(Diagram):
    """
    A plot of data with 2D spatial information, either world coordinates or instrument
    coordinates that are connected with a projection.
    """

    def __init__(self, axes, projection=None, orig=None):
        super().__init__(axes, orig=orig)

        if not isinstance(orig, SpatialDiagram):
            self.__projection = projection
        else:
            self.__projection = projection or orig.__projection

    def __get_projection(self):
        return self.__projection

    projection = property(__get_projection)

    def apply(self, ax):
        super().apply(ax)

    def project_coords(self, ax: plt.Axes, *coords, native_frame=None):
        ctype, coords = normalize_coords(*coords)
        need_conversion = native_frame != self._get_native_frame()

        if need_conversion:
            if native_frame == 'pixel':
                coords, _ = self.projection.pixel_to_world(coords)
            elif  native_frame == 'world':
                coords, _ = self.projection.world_to_pixel(coords)
            else:
                raise NotImplementedError()

        return denormalize_coords(ctype, coords)

    def _get_coords(self, ax: plt.Axes, *coords, mask_fov=False, mask=None, s=None, native_frame=None, **kwargs):
        """
        Compare data frames (world or pixel), and perform conversion, if necessary.
        Also set the transformation of the axes, if WCS projection is used.
        """

        ctype, coords = normalize_coords(*coords)
        need_conversion = native_frame != self._get_native_frame()

        if mask_fov or need_conversion:
            if native_frame == 'pixel':
                cc, fov_mask = self.projection.pixel_to_world(coords, mask=mask)
            elif  native_frame == 'world':
                cc, fov_mask = self.projection.world_to_pixel(coords, mask=mask)
            else:
                raise NotImplementedError()

            if need_conversion:
                coords = cc
        else:
            fov_mask = None

        if isinstance(ax, WCSAxes):
            transform = kwargs.pop('transform', ax.get_transform(native_frame))
        else:
            transform = None

        return ctype, coords, mask, s, fov_mask, transform

    def _get_native_frame(self, native_frame=None):
        raise NotImplementedError()

    def plot(self, ax: plt.Axes, *coords, fmt=None, native_frame=None,
             mask_fov=False, mask=None, s=None, **kwargs):

        ctype, coords, mask, s, fov_mask, transform = self._get_coords(ax, 
            *coords, mask_fov=mask_fov, mask=mask, s=s, native_frame=native_frame, **kwargs)
    
        l = super().plot(ax, coords[..., 0], coords[..., 1], fmt=fmt, transform=transform, mask=mask, s=s, **kwargs)
        return l

    def scatter(self, ax: plt.Axes, *coords, native_frame=None,
                mask_fov=True, mask=None, s=None, **kwargs):

        ctype, coords, mask, s, fov_mask, transform = self._get_coords(ax, 
            *coords, mask_fov=mask_fov, mask=mask, s=s, native_frame=native_frame, **kwargs)
        
        l = super().scatter(ax, coords[..., 0], coords[..., 1], transform=transform, s=s, **kwargs)
        return l

    def plot_catalog(self, ax: plt.Axes, catalog, population_id=None, apply_categories=None, g=None,
                     observed=None, mask_fov=False, mask=None, s=None, **kwargs):
        
        if isinstance(catalog, Observation):
            l = self.plot_observation(ax, catalog, observed=observed, mask_fov=mask_fov, mask=mask, s=s, **kwargs)
        elif isinstance(catalog, Simulation):
            observed = observed if observed is not None else True
            apply_categories = apply_categories if apply_categories is not None else True

            l = self.plot_simulation(ax, catalog, population_id=population_id, apply_categories=apply_categories, g=g,
                observed=observed, mask_fov=mask_fov, mask=mask, s=s, **kwargs)
        else:
            raise NotImplementedError()

        self.apply(ax)

        return l
        

    def plot_simulation(self, ax: plt.Axes, catalog, population_id=None, apply_categories=False, g=None,
            observed=None, mask_fov=False, mask=None, s=None, **kwargs):
        
        # TODO: implement
        raise NotImplementedError()

    def plot_observation(self, ax: plt.Axes, catalog, 
            observed=None, mask_fov=False, mask=None, s=None, **kwargs):

        observed = observed if observed is not None else catalog.observed
        ra, dec = catalog.get_coords(mask=mask)

        # Do not pass on mask
        return self.scatter(ax, ra, dec, mask_fov=mask_fov, mask=False, s=s, native_frame='world', **kwargs)