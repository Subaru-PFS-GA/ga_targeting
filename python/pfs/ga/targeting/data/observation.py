import logging
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from ..util import *
from ..photometry import Photometry, Color, Magnitude
from .catalog import Catalog

class Observation(Catalog):
    def __init__(self, data: pd.DataFrame = None, name=None, frame='icrs', equinox='J2000', orig=None):
        super().__init__(name=name, frame=frame, equinox=equinox, orig=orig)

        if not isinstance(orig, Observation):
            self.__data: pd.DataFrame = data
        else:
            self.__data: pd.DataFrame = data if data is not None else safe_deep_copy(orig.__data)

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

    def get_data(self, mask=None, func=None, selection=None) -> pd.DataFrame:
        """
        Return all or a subset of the catalog as a data frame. Only one of the
        arguments `mask`, `func` and `selection` can be specified at a time.

        Arguments
        ---------
        mask : array-like
            Boolean mask to select rows.
        func : callable
            Function to apply to the data to select rows.
        selection : Selection
            Selection object to apply to the data to select rows.
        """

        not_none = 0
        for i in [mask, func, selection]:
            if i is not None:
                not_none += 1

        if not_none > 1:
            raise RuntimeError('Only one of `mask`, `func` and `selection` can be not None.')
        elif mask is not None:
            return self.__data[mask]
        elif func is not None:
            return func(self.__data)
        elif selection is not None:
            mask = selection.apply(self)
            return self.__data[mask]
        else:
            return self.__data
        
    def filter(self, *, mask=None, func=None, selection=None):
        """
        Return all or a subser of the catalog as a new catalog. Only one of the
        arguments `mask`, `func` and `selection` can be specified at a time.

        Arguments
        ---------
        mask : array-like
            Boolean mask to select rows.
        func : callable
            Function to apply to the data to select rows.
        selection : Selection
            Selection object to apply to the data to select rows.
        """

        data = self.get_data(mask=mask, func=func, selection=selection)
        return type(self)(data=data, orig=self)

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

    def calculate_flux(self, unit=u.nJy, force=False):
        """
        For each magnitude in the catalog, calculate the flux and its error, if
        not already available.

        If `force` is True, then the flux is recalculated even if it is already available.

        Parameters
        ----------
        force : bool
            If True, recalculate the flux even if it is already available.
        """

        for p in self.photometry.values():
            for m in p.magnitudes.values():
                if m.get_name('obs_') in self.__data and (force or m.get_name('obs_flux_') not in self.__data):
                    mag, mag_err = self.get_magnitude(m, observed=True, dered=False)
                    
                    if mag is not None:
                        # Convert to flux in nJy using astropy
                        flux = (np.array(mag) * u.ABmag).to_value(unit)
                        self.__data[m.get_name('obs_flux_')] = flux

                        if mag_err is not None:
                            flux_err = flux * np.array(mag_err) * np.log(10) / 2.5
                            self.__data[m.get_name('err_flux_')] = flux_err
