import logging
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

from ..util.args import *
from ..util.coords import *
from ..util import safe_deep_copy, ReadOnlyDict
from ..photometry import Photometry, Color, Magnitude
from ..diagram import ColorAxis, MagnitudeAxis, DiagramValueProvider

from ..setup_logger import logger

class Catalog(DiagramValueProvider):
    def __init__(self, name=None, frame='icrs', equinox='J2000', orig=None):
        if not isinstance(orig, Catalog):
            self.__name = name
            self.__photometry = {}
            self.__frame = frame
            self.__equinox = equinox
        else:
            self.__name = name or orig.__name
            self.__photometry = safe_deep_copy(orig.__photometry)
            self.__frame = frame or orig.__frame
            self.__equinox = equinox or orig.__equinox

    def __len__(self):
        # TODO: return number of simulated stars per population
        raise NotImplementedError()

    def __get_name(self) -> str:
        return self.__name
    
    def __set_name(self, value: str):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_photometry(self):
        return ReadOnlyDict(self.__photometry)
    
    def _set_photometry(self, value):
        self.__photometry = value
    
    photometry = property(__get_photometry)

    def __get_frame(self):
        return self.__frame

    def __set_frame(self, value):
        self.__frame = value

    frame = property(__get_frame, __set_frame)

    def __get_equinox(self):
        return self.__equinox
    
    def __set_equinox(self, value):
        self.__equinox = value

    equinox = property(__get_equinox, __set_equinox)

    def __get_observed(self):
        raise NotImplementedError()

    observed = property(__get_observed)

    def append_photometry(self, photometry):
        self.__photometry[photometry.name] = photometry

    def get_magnitude(self, magnitude: Magnitude, observed=False, mask=None):
        raise NotImplementedError()

    def get_color(self, color: Color, observed=False, mask=None):
        m1, s_m1 = self.get_magnitude(color.magnitudes[0], observed=observed, mask=mask)
        m2, s_m2 = self.get_magnitude(color.magnitudes[1], observed=observed, mask=mask)

        # TODO: correlated error?
        if s_m1 is not None and s_m2 is not None:
            s_c = np.sqrt(s_m1**2 + s_m2**2)
        else:
            s_c = None

        return m1 - m2, s_c

    def get_id(self, id=None, mask=None):
        mask = np.s_[()] if mask is None else mask

        id = np.array(self.data[id or 'objid'])[mask]
        return id
    
    def has_coords(self, ra=None, dec=None):
        ra = ra if ra is not None else 'RA'
        dec = dec if dec is not None else 'Dec'

        return ra in self.data and dec in self.data

    def get_coords(self, ra=None, dec=None, mask=None, ctype=None):
        ra = ra if ra is not None else 'RA'
        dec = dec if dec is not None else 'Dec'
        
        mask = np.s_[()] if mask is None else mask

        ra = np.array(self.data[ra])[mask]
        dec = np.array(self.data[dec])[mask]

        if ctype is not None:
            _, coords = normalize_coords(ra, dec)
            return denormalize_coords(ctype, coords)
        else:
            return ra, dec

    def get_skycoords(self, mask=None, frame=None, equinox=None, **kwargs):

        # TODO: what if we have errors?

        frame = frame or self.__frame
        equinox = equinox or self.__equinox

        ra, dec = self.get_coords(mask=mask)
        coords = SkyCoord(ra, dec, unit=u.degree, frame=frame, equinox=equinox, **kwargs)
        return coords
    
    def cone_search(self, pos, rad):
        pos = normalize_pos(pos)
        rad = normalize_angle(rad, u.arcmin)

        coords = self.get_skycoords()
        d = pos.separation(coords)
        mask = (d.arcmin <= rad.arcmin)

        return mask
    
    def random_sample(self, p):
        n = len(self)
        c = np.random.choice(np.arange(n), int(n * p), replace=False)
        mask = np.full((n,), False)
        mask[c] = True

        return mask
    
    def cross_match(self, other, max_separation=1, mask=None, mask_other=None):
        """
        Cross match two catalogs.
        
        Use built-in astropy functionality to perform the matching.
        """

        max_separation = normalize_angle(max_separation, u.arcsec, allow_none=True)

        c1 = self.get_skycoords(mask=mask)
        c2 = other.get_skycoords(mask=mask_other)

        idx, separation, _ = c1.match_to_catalog_sky(c2)
        
        # Apply cut on separation
        if max_separation is not None:
            idx[separation > max_separation] = -1

        return idx, separation