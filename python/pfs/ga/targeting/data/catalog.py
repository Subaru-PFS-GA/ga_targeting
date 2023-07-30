import logging
import numpy as np
import pandas as pd

from ..util.coords import *
from ..util import safe_deep_copy, ReadOnlyDict
from ..photometry import Photometry, Color, Magnitude
from ..diagram import ColorAxis, MagnitudeAxis, DiagramValueProvider

class Catalog(DiagramValueProvider):
    def __init__(self, name=None, orig=None):
        if not isinstance(orig, Catalog):
            self.__name = name
            self.__photometry = {}
        else:
            self.__name = name or orig.__name
            self.__photometry = safe_deep_copy(orig.__photometry)

    def __len__(self):
        # TODO: return number of simulated stars per population
        raise NotImplementedError()

    def __get_name(self) -> str:
        # TODO: return number of simulated times populations
        return self.__name
    
    def __set_name(self, value: str):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_photometry(self):
        return ReadOnlyDict(self.__photometry)

    photometry = property(__get_photometry)

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

    def get_coords(self, ra=None, dec=None, mask=None, ctype=None):
        mask = np.s_[()] if mask is None else mask

        ra = np.array(self.data[ra or 'RA'])[mask]
        dec = np.array(self.data[dec or 'Dec'])[mask]

        if ctype is not None:
            _, coords = normalize_coords(ra, dec)
            return denormalize_coords(ctype, coords)
        else:
            return ra, dec

        