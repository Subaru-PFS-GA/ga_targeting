import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.coordinates import Angle

from ..util import *

from .projection import Projection

class WcsProjection(Projection):
    def __init__(self, pointing=None, scale=None, proj=None, orig=None):
        super().__init__(pointing=pointing, orig=orig)

        if not isinstance(orig, WcsProjection):
            self.__scale = scale or 1.0
            self.__proj = proj or 'TAN'
        else:
            self.__scale = scale or orig.__scale
            self.__proj = proj or orig.__proj

    def __get_scale(self):
        return self.__scale
    
    def __set_scale(self, value):
        self.__scale = value

    scale = property(__get_scale, __set_scale)

    def __get_proj(self):
        return self.__proj
    
    def __set_proj(self, value):
        self.__proj = value

    proj = property(__get_proj, __set_proj)

    def __get_wcs(self):
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["RA---{}".format(self.__proj), "DEC--{}".format(self.__proj)]
        w.wcs.crpix = [0, 0]
        w.wcs.crval = [self.pointing.ra, self.pointing.dec]
        w.wcs.cdelt = [self.__scale, self.__scale]
        w.wcs.set_pv([(2, 1, 0.0)])     # TODO: set position angle

        return w

    wcs = property(__get_wcs)

    def world_to_pixel(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)
        ra, dec = split_coords(coords)
        w = self.__get_wcs()
        x, y = w.wcs_world2pix(ra, dec, 1)
        fov_mask = np.full_like(x, True, dtype=bool)
        return denormalize_coords(ctype, x, y), fov_mask

    def pixel_to_world(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)
        w = self.__get_wcs()
        x, y = split_coords(coords)
        ra, dec = w.wcs_pix2world(x, y, 1)
        fov_mask = np.full_like(ra, True, dtype=bool)
        return denormalize_coords(ctype, stack_coords(ra, dec)), fov_mask
