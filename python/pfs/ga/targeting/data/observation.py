import logging
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from ..util import *
from ..photometry import Photometry, Color, Magnitude
from .catalog import Catalog

class Observation(Catalog):
    def __init__(self, name=None, frame='icrs', equinox='J2000', orig=None):
        super().__init__(name=name, frame=frame, equinox=equinox, orig=orig)

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

    def get_data(self, mask=None, filter=None):
        if mask is not None and filter is not None:
            raise RuntimeError('Either `mask` or `filter` can be not None.')
        elif mask is not None:
            return self.__data[mask]
        elif filter is not None:
            return filter(self.__data)
        else:
            return self.__data

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

    def append_columns(self, other, other_key, source_columns, target_columns=None):
        """
        Merge columns `source_columns` from `other` by joining the column of self 'objid'
        on the column ˙other_key` of `other`. The key `other_key` is not copied.
        """
        this_columns = [ 'objid' ]
        other_columns = [ other_key ] + source_columns
        df = self.__data[this_columns].join(other[other_columns].set_index(other_key), on='objid', how='left')
        for s, t in zip(source_columns, target_columns if target_columns is not None else source_columns):
            self.__data[t] = df[s]