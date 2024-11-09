import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from ..setup_logger import logger
from ..util import *
from ..util.pandas import *
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

    def __get_columns(self):
        return self.__data.columns
    
    columns = property(__get_columns)

    def get_data(self, mask=None, filter=None, selection=None) -> pd.DataFrame:
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
        for i in [mask, filter, selection]:
            if i is not None:
                not_none += 1

        if not_none > 1:
            raise RuntimeError('Only one of `mask`, `func` and `selection` can be not None.')
        elif mask is not None:
            return self.__data[mask]
        elif filter is not None:
            return filter(self.__data)
        elif selection is not None:
            mask = selection.apply(self)
            return self.__data[mask]
        else:
            return self.__data
        
    def filter(self, *, mask=None, filter=None, selection=None):
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

        data = self.get_data(mask=mask, filter=filter, selection=selection)
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
                        logger.warning(f'Extinction correction is not available in catalog `{self.name}` for magnitude `{magnitude.filter}`.')
            else:
                mag = get_value(magnitude.get_name('obs_'), mask)
        else:
            if dered:
                mag = get_value(magnitude.get_name(), mask)
                ext = get_value(magnitude.get_name('ext_'), mask)
                if ext is not None:
                    mag = mag - ext
                else:
                    logger.warning(f'Extinction correction is not available in catalog `{self.name}` for magnitude `{magnitude.filter}`.')
            else:
                mag = get_value(magnitude.get_name(), mask)

        if mag is None:
            raise Exception(f'Magnitude `{magnitude.filter}` not available for the catalog `{self.name}`')

        err = get_value(magnitude.get_name('err_'), mask)

        return mag, err
    
    def append_column(self, name, values, dtype):
        """
        Append a column to the catalog.
        """
        
        pd_append_column(self.__data, name, values, dtype=dtype)

    def set_column_dtype(self, name, dtype):
        """
        Set the data type of a column.
        """
        
        self.__data[name] = self.__data[name].astype(dtype)

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

    def merge(self, other, idx, columns, mask=None, mask_other=None, other_prefix=None):
        """
        Merge two catalogs based on indices pointing to the other. This will
        result in a left outer join type match.
        """

        ix1 = self.data.index[mask]
        ix2 = other.data.index[idx[mask]]
        
        for c in columns:
            self.data[c] = np.nan
            self.data.loc[ix1, c] = np.array(other.data.loc[ix2, c])
                
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

    def calculate_flux_filters(self, filters, unit=u.nJy):
        """
        Calculate flux in filters. This is necessary to fill in missing data for
        the design files.
        """

        # TODO: column names are inconsistent with the version of this function that
        #       works with photometry objects

        for filter, columns in filters.items():
            any_flux = None         # First available flux value for a filter
            any_flux_err = None     # First available flux error for a filter
            any_flux_key = None     # Key of the first available flux column
            for prefix in ['', 'psf_', 'fiber_', 'total_']:
                flux_key = prefix + 'flux'
                flux_err_key = prefix + 'flux_err'
                mag_key = prefix + 'mag'
                mag_err_key = prefix + 'mag_err'

                flux_col_canonical = f'{filter}_{prefix}flux'
                flux_err_col_canonical = f'{filter}_{prefix}flux_err'

                flux_col = columns[flux_key] if flux_key in columns else None
                flux_err_col = columns[flux_err_key] if flux_err_key in columns else None
                mag_col = columns[mag_key] if mag_key in columns else None
                mag_err_col = columns[mag_err_key] if mag_err_key in columns else None

                # Calculate flux from the magnitude
                flux = None
                flux_err = None
                if flux_col is None:
                    logger.warning(f'Missing flux column `{flux_key}` in target list {self.name}.')

                    if mag_col is not None:
                        # The magnitude is available to calculate the flux from
                        logger.info(f'Calculating `{flux_key}` from `{mag_key}` in filter {filter}.')

                        mag = self.__data[mag_col]
                        mag_err = self.__data[mag_err_col] if mag_err_col is not None else None

                        flux = (np.array(mag) * u.ABmag).to_value(unit)
                        if mag_err is not None:
                            flux_err = 0.4 * np.log(10) * flux * np.array(mag_err)
                    elif any_flux is not None:
                        # No magnitude is available, copy flux with another prefix
                        logger.info(f'Copying `{flux_key}` from `{any_flux_key}` in filter {filter}.')
                        flux = any_flux
                        flux_err = any_flux_err
                else:
                    flux = self.__data[flux_col]
                    if flux_err_col is not None:
                        flux_err = self.__data[flux_err_col]

                # Save the flux to the target list with canonical column name
                self.__data[flux_col_canonical] = flux if flux is not None else np.nan
                self.__data[flux_err_col_canonical] = flux_err if flux_err is not None else np.nan

                if any_flux is None:
                    any_flux = self.__data[flux_col_canonical]
                    any_flux_err = self.__data[flux_err_col_canonical]
                    any_flux_key = flux_key
    
    def calculate_flux_bands(self, bands, unit=u.nJy):
        """
        Calculate flux in filters. This is necessary to fill in missing data for
        the design files.
        """

        for band, columns in bands.items():
            any_flux = None
            any_flux_err = None
            any_flux_key = None
            for prefix in ['', 'psf_', 'fiber_', 'total_']:
                # Column keys
                flux_key = prefix + 'flux'
                flux_err_key = prefix + 'flux_err'
                mag_key = prefix + 'mag'
                mag_err_key = prefix + 'mag_err'

                # Column names
                flux_col = columns[flux_key] if flux_key in columns else None
                flux_err_col = columns[flux_err_key] if flux_err_key in columns else None
                mag_col = columns[mag_key] if mag_key in columns else None
                mag_err_col = columns[mag_err_key] if mag_err_key in columns else None

                flux_col_canonical = f'{prefix}flux_{band}'
                flux_err_col_canonical = f'{prefix}flux_err_{band}'

                # Calculate flux from the magnitude
                flux = None
                flux_err = None
                if flux_col is None:
                    logger.warning(f'Missing flux column `{flux_key}` in target list {self.name}.')

                    if mag_col is not None:
                        # The magnitude is available to calculate the flux from
                        logger.info(f'Calculating `{flux_key}` from `{mag_key}` in filter {filter}.')

                        mag = self.__data[mag_col]
                        mag_err = self.__data[mag_err_col] if mag_err_col is not None else None

                        flux = (np.array(mag) * u.ABmag).to_value(unit)
                        if mag_err is not None:
                            flux_err = 0.4 * np.log(10) * flux * np.array(mag_err)
                    elif any_flux is not None:
                        # No magnitude is available, copy flux with another prefix
                        logger.info(f'Copying `{flux_key}` from `{any_flux_key}` in filter {filter}.')
                        flux = any_flux
                        flux_err = any_flux_err
                else:
                    flux = self.__data[flux_col]
                    if flux_err_col is not None:
                        flux_err = self.__data[flux_err_col]

                # Save the flux to the target list with canonical column name
                self.__data[flux_col_canonical] = flux if flux is not None else np.nan
                self.__data[flux_err_col_canonical] = flux_err if flux_err is not None else np.nan

                if any_flux is None:
                    any_flux = self.__data[flux_col_canonical]
                    any_flux_err = self.__data[flux_err_col_canonical]
                    any_flux_key = flux_key