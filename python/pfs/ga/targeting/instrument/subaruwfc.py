import numpy as np
from collections import namedtuple
from matplotlib.patches import Circle

from ics.cobraOps.Bench import Bench
from pfs.utils.coordinates.CoordTransp import CoordinateTransform

from pfs.ga.common.util import *
from pfs.ga.common.diagram import styles

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
    

    def world_to_fp_pos(self, ra, dec, pmra=None, pmdec=None, parallax=None, epoch=2015.5, mask=None):
        # TODO: use mask
        # mask = mask if mask is not None else np.full(ra.shape[:-1], True, dtype=bool)

        # The center of input coordinates in the same unit as xyin.
        cent = np.array([[ self.pointing.ra ], [ self.pointing.dec ]])
                        
        # Position angle in unit of degree for sky_pfi* transformation.
        pa = self.pointing.posang

        # Input coordinates. Namely. (Ra, Dec) in unit of degree for sky with shape (2, N)
        _, coords = normalize_coords(ra, dec)
        xyin = coords.T

        # The proper motion of the targets used for sky_pfi transformation.
        # The unit is mas/yr, the shape is (2, N)
        if pmra is not None and pmdec is not None:
            pm = np.stack([ pmra, pmdec ], axis=0)
        else:
            pm = None

        # The parallax of the cordinatess used for sky_pfi transformation.
        # The unit is mas, the shape is (1, N)
        # if parallax is not None:
        #     parallax = parallax[None, :]
        # else:
        #     parallax = None

        # Observation time UTC in format of %Y-%m-%d %H:%M:%S
        obs_time = self.pointing.obs_time.to_value('iso')

        fp_pos = CoordinateTransform(xyin=xyin,
                                     mode="sky_pfi",
                                     # za=0.0, inr=0.0,     # These are overriden by function
                                     cent=cent,
                                     pa=pa,
                                     pm=pm,
                                     par=parallax,
                                     time=obs_time,
                                     epoch=epoch)
                
        xy = fp_pos[:2, :].T
        
        # Old version of sky_pfi converted returned the radius, now we have to calculate it
        # r = fp_pos[3, :]
        r = np.sqrt(np.sum(xy ** 2, axis=-1))

        # Mask coordinates inside fov radius
        fov_mask = (r <= (498 / 2.0))
    
        return xy, fov_mask

    def world_to_pixel(self, *coords, mask=None):
        # This is a simplified signature to be used for plotting only

        ctype, coords = normalize_coords(*coords)

        xy, fov_mask = self.world_to_fp_pos(coords[..., 0], coords[..., 1], mask=mask)
    
        return denormalize_coords(ctype, xy), fov_mask


    def pixel_to_world(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)

        # Make sure coordinates are 2d
        if coords.ndim < 2:
            coords = coords[None, :]

        # TODO: use mask
        mask = mask if mask is not None else np.full_like(coords.T, True, dtype=bool)

        # Different versions of pfs_utils require different shape for
        # cent.

        sky_pos = CoordinateTransform(xyin=coords.T, mode="pfi_sky",
            za=0.0, inr=0.0,
            cent=np.array([[self.pointing.ra, self.pointing.dec]]).T,
            # cent=np.array([self.pointing.ra, self.pointing.dec]),
            pa=self.pointing.posang,
            time=self.pointing.obs_time.to_string())

        radec = sky_pos[:2, :].T
        r = coords[..., 0]**2 + coords[..., 1]**2
        fov_mask = (r <= (498 / 2.0)**2)

        return denormalize_coords(ctype, radec), fov_mask

    def __get_outline_points_pixel(self, res=None):
        res = res if res is not None else 360

        phi = np.linspace(0, 2 * np.pi, res)
        xy = 498 / 2.0 * np.stack([np.cos(phi), np.sin(phi)], axis=-1)

        return xy, np.full(phi.shape, True, dtype=bool)

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