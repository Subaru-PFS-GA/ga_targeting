import numpy as np
from collections import namedtuple
from matplotlib.patches import Circle

from ics.cobraOps.Bench import Bench
from pfs.utils.coordinates.CoordTransp import CoordinateTransform

from ..util import *
from ..diagram import styles
from ..projection.telescopeprojection import TelescopeProjection
from . instrument import Instrument

class SubaruWFC(Instrument, TelescopeProjection):
    """
    Focal plane projection of the Subaru telescope with the Wide Field Corrector.
    """

    def __init__(self, pointing=None, orig=None):
        Instrument.__init__(self, orig=orig)
        TelescopeProjection.__init__(self, pointing=pointing, orig=orig)

        if not isinstance(orig, SubaruWFC):
            pass
        else:
            pass

    def world_to_pixel(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)

        # TODO: use mask
        mask = mask if mask is not None else np.full(coords.shape[:-1], True, dtype=bool)

        fp_pos = CoordinateTransform(xyin=coords.T, mode="sky_pfi",
            za=0.0, inr=0.0,
            cent=np.array([[ self.pointing.ra ], [ self.pointing.dec ]]),
            pa=self.pointing.posang,
            time=self.pointing.time)
        
        xy = fp_pos[:2, :].T
        
        # Old version of sky_pfi converted returned the radius, now we have to calculate it
        # r = fp_pos[3, :]
        r = np.sqrt(np.sum(xy ** 2, axis=-1))

        # Mask coordinates inside fov radius
        fov_mask = (r <= (498 / 2.0))
    
        return denormalize_coords(ctype, xy), fov_mask

    def pixel_to_world(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)

        # TODO: use mask
        mask = mask if mask is not None else np.full_like(coords.T, True, dtype=bool)

        sky_pos = CoordinateTransform(xyin=coords.T, mode="pfi_sky",
            za=0.0, inr=0.0,
            cent=np.array([self.pointing.ra, self.pointing.dec]), pa=self.pointing.posang,
            time=self.pointing.time)

        radec = sky_pos[:2, :].T
        r = coords[..., 0]**2 + coords[..., 1]**2
        fov_mask = (r <= (498 / 2.0)**2)

        return denormalize_coords(ctype, radec), fov_mask

    def __get_outline_points_pixel(self, res=None):
        res = res if res is not None else 360

        phi = np.linspace(0, 2 * np.pi, res)
        xy = 498 / 2.0 * np.stack([np.cos(phi), np.sin(phi)], axis=-1)

        return xy, np.full(phi.shape, True, dtype=np.bool)

    def __get_outline_points_world(self, res=None):
        xy, mask = self.__get_outline_points_pixel(res=res)
        radec, mask = self.pixel_to_world(xy, mask=mask)
        return radec, mask

    def plot_focal_plane(self, ax, fp, res=None, **kwargs):
        style = styles.solid_line(**kwargs)
        
        xy, mask = self.__get_outline_points_pixel(res=res)
        fp.plot(ax, xy[..., 0], xy[..., 1], native_frame='pixel', **style)

    def plot_field_of_view(self, ax, fov, res=None, **kwargs):
        style = styles.solid_line(**kwargs)

        xy, mask = self.__get_outline_points_world(res=res)
        fov.plot(ax, xy[..., 0], xy[..., 1], mask_fov=False, mask=mask, native_frame='world', **style)