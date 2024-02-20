import logging
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from ..util import *
from ..photometry import Photometry, Color, Magnitude
from .catalog import Catalog

class Observation(Catalog):
    def __init__(self, name=None, orig=None):
        super(Observation, self).__init__(name=name, orig=orig)

        if not isinstance(orig, Observation):
            self.__data: pd.DataFrame = None
        else:
            self.__data = safe_deep_copy(orig.__data)

    def __len__(self):
        return len(self.__data)

    def __get_shape(self):
        return (len(self.__data),)

    shape = property(__get_shape)

    def __get_data(self) -> pd.DataFrame:
        return self.__data

    def _set_data(self, data: pd.DataFrame):
        self.__data = data

    data = property(__get_data)

    def __get_observed(self):
        return True

    observed = property(__get_observed)

    def has_magnitude(self, magnitude: Magnitude, observed=True, dered=True):
        if observed:
            if dered:
                if magnitude.get_name('dered_') in self.__data:
                    return True
                if magnitude.get_name('obs_') in self.__data and magnitude.get_name('ext_') in self.__data:
                    return True
            else:
                if magnitude.get_name('obs_') in self.__data:
                    return True
        else:
            if dered:
                if magnitude.get_name() in self.__data and magnitude.get_name('ext_'):
                    return True
            else:
                if magnitude.get_name() in self.__data:
                    return True
                    
        return False

    def get_magnitude(self, magnitude: Magnitude, observed=True, dered=True, mask=None):
        mask = mask if mask is not None else slice(None)

        def get_value(k, mask):
            if k in self.__data:
                return self.__data[k][mask]
            else:
                return None

        mag = None
        ext = None
        err = None

        if observed:
            if dered:
                mag = get_value(magnitude.get_name('dered_'), mask)
            
                if mag is None:
                    mag = get_value(magnitude.get_name('obs_'), mask)
                    ext = get_value(magnitude.get_name('ext_'), mask)
                    if ext is not None:
                        mag = mag - ext
                    else:
                        logging.warning(f'Extinction correction is not available in catalog `{self.name}` for magnitude `{magnitude.filter}`.')
            else:
                mag = get_value(magnitude.get_name('obs_'), mask)
        else:
            if dered:
                mag = get_value(magnitude.get_name(), mask)
                ext = get_value(magnitude.get_name('ext_'), mask)
                if ext is not None:
                    mag = mag - ext
                else:
                    logging.warning(f'Extinction correction is not available in catalog `{self.name}` for magnitude `{magnitude.filter}`.')
            else:
                mag = get_value(magnitude.get_name(), mask)

        if mag is None:
            raise Exception(f'Magnitude `{magnitude.filter}` not available for the catalog `{self.name}`')

        err = get_value(magnitude.get_name('err_'), mask)

        return mag, err

    def cross_match(self, other, mask=None, mask_other=None):
        ra1, dec1 = self.get_coords(mask=mask)
        c1 = SkyCoord(ra1 * u.deg, dec1 *u.deg)

        ra2, dec2 = other.get_coords(mask=mask_other)
        c2 = SkyCoord(ra2 * u.deg, dec2 *u.deg)

        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        
        # idx: matching ids in `other`
        # d2d: distance 2d (angle)
        # d3d: distance 3d (?)

        return idx, d2d.value