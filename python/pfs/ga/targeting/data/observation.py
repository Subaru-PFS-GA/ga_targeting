import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

from ..setup_logger import logger
from ..util import *
from ..util.pandas import *
from ..photometry import Photometry, Color, Magnitude
from .catalog import Catalog

class Observation(Catalog):
    def __init__(self, data: pd.DataFrame = None, name=None, frame='icrs', equinox=None, epoch=None, orig=None):
        super().__init__(name=name, frame=frame, equinox=equinox, epoch=epoch, orig=orig)

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
        filter : callable
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
            return self.__data.loc[mask]
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
        data.reset_index(inplace=True, drop=True)
        return type(self)(data=data, orig=self)

    def apply_magnitude_limits(self, limits):

        def get_magnitude(key):
            # Split the filter name into magnitude_type, photometric_system, filter_name
            parts = key.split('_')
            if len(parts) > 1:
                [p, m] = parts[-2:]
                magnitude_type = parts[:-2]
                mag, mag_err = self.get_magnitude(self.photometry[p].magnitudes[m], magnitude_type=magnitude_type)
                return mag, mag_err
            else:
                raise ValueError(f"Invalid filter: {key}")

        def append_mask(mask, new_mask):
            if mask is None:
                return new_mask
            else:
                return mask & new_mask

        mask = None
        for key, values in limits.items():
            # Split the key name at dashes to see if it is a color definition
            parts = key.split('-')
            if len(parts) == 1:
                mag, _ = get_magnitude(parts[0])
                if mag is not None:
                    mask = append_mask(mask, (mag >= values[0]) & (mag <= values[1]))
                else:
                    raise ValueError(f"Invalid filter: {parts[0]}")
            elif len(parts) == 2:
                mag1, _ = get_magnitude(parts[0].strip())
                mag2, _ = get_magnitude(parts[1].strip())
                if mag1 is not None and mag2 is not None:
                    mask = append_mask(mask, (mag1 - mag2 >= values[0]) & (mag1 - mag2 <= values[1]))
                else:
                    raise ValueError(f"Invalid filter: {parts[0]} or {parts[1]}")
            else:
                raise ValueError(f"Invalid filter or color definition: {key}")

        if mask is not None:
            return self.filter(mask=mask)

    def __get_observed(self):
        return True

    observed = property(__get_observed)

    def __match_magnitude(self, magnitude: Magnitude):
        p = magnitude.photometry.name
        m = magnitude.filter
        if p in self.photometry and m in self.photometry[p].magnitudes:
            return self.photometry[p].magnitudes[m]
        else:
            return magnitude

    def has_magnitude(self, magnitude: Magnitude, observed=True, dered=True):
        magnitude = self.__match_magnitude(magnitude)

        if magnitude.columns is not None and len(magnitude.columns) > 0:
            return self.has_magnitude_by_column_name(magnitude, observed=observed, dered=dered)
        else:
            return self.has_magnitude_by_name(magnitude, observed=observed, dered=dered)

    def has_magnitude_by_name(self, magnitude: Magnitude, observed=True, dered=True):
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
    
    def has_magnitude_by_column_name(self, magnitude: Magnitude, observed=True, dered=True):
        """
        Return True if the magnitude in the specified filter is available in any form or
        can be calculated from a corresponding flux value.
        """

        mag, mag_err, flux, flux_err, ext, magnitude_type = self.get_magnitude_column_names(magnitude)
        return mag is not None or flux is not None
    
    def get_magnitude(self, magnitude: Magnitude, observed=None, dered=None, magnitude_type=None, mask=None):
        magnitude = self.__match_magnitude(magnitude)

        if magnitude.columns is not None and len(magnitude.columns) > 0:
            return self.get_magnitude_by_column(magnitude, observed=observed, dered=dered, magnitude_type=magnitude_type, mask=mask)
        else:
            return self.get_magnitude_by_name(magnitude, observed=observed, dered=dered, magnitude_type=magnitude_type, mask=mask)
        
    def get_magnitude_column_names(self, magnitude: Magnitude, magnitude_type=None):
        """
        Return the data frame columns corresponding to the magnitude, flux and their errors.
        The type of magnitude is selected in the precendence order defined in `magnitude_type`.
        """

        # Normalize the magnitude type precedence order parameter
        if magnitude_type is not None and not isinstance(magnitude_type, Iterable):
            magnitude_type = [ magnitude_type ]
        elif magnitude_type is None \
             or isinstance(magnitude_type, Iterable) and len(magnitude_type) == 0:
            
            magnitude_type = [ 'psf', 'fiber', 'total', '' ]

        # Try to find the requested magnitude columns by descreasing order of precedence
        for t in magnitude_type:
            prefix = f'{t}_' if t != '' else ''

            mag = None
            mag_err = None
            flux = None
            flux_err = None
            ext = None

            # Get the column names for the magnitude and its error
            # Alternatively, use the flux column and its error if the magnitude is not available
            mag = magnitude.columns.get(f'{prefix}mag')
            mag_err = magnitude.columns.get(f'{prefix}mag_err')
            flux = magnitude.columns.get(f'{prefix}flux')
            flux_err = magnitude.columns.get(f'{prefix}flux_err')
            ext = magnitude.columns.get(f'{prefix}mag_ext')
            
            # Reset values that are not available in the data frame
            def reset(k):
                return k if k is not None and k in self.__data else None
            
            mag = reset(mag)
            mag_err = reset(mag_err)
            flux = reset(flux)
            flux_err = reset(flux_err)
            ext = reset(ext)

            # Error columns cannot exist without the value column
            if mag is None:
                mag_err = None
                ext = None

            if flux is None:
                flux_err = None

            if mag is not None or flux is not None:
                return mag, mag_err, flux, flux_err, ext, magnitude_type

        return None, None, None, None, None, None
            
    def get_magnitude_by_column(self, magnitude: Magnitude, observed=None, dered=None, magnitude_type=None, mask=None):
        """
        Return the values of a magnitude from the data frame. If the magnitude is
        not available, but there is a flux column, do the conversion on-the-fly.
        """

        observed = observed if observed is not None else True
        dered = dered if dered is not None else True
        mask = mask if mask is not None else slice(None)

        if isinstance(magnitude_type, str):
            magnitude_type = [ magnitude_type ]
        elif magnitude_type is not None and not isinstance(magnitude_type, Iterable):
            magnitude_type = [ magnitude_type ]
        elif magnitude_type is None \
             or isinstance(magnitude_type, Iterable) and len(magnitude_type) == 0:
            
            magnitude_type = [ 'psf', 'fiber', 'total', '' ]

        def get_value(k, mask):
            if k is not None and k in self.__data:
                return self.__data[k][mask]
            else:
                return None
            
        # Depending on what's available in the column list, return the magnitude
        # Ignore the 'observed' flag
        mag, mag_err, flux, flux_err, ext, magnitude_type = \
            self.get_magnitude_column_names(magnitude, magnitude_type=magnitude_type)
        
        # Get the columns from the data frame
        mag = get_value(mag, mask)
        mag_err = get_value(mag_err, mask)
        ext = get_value(ext, mask)
        flux = get_value(flux, mask)
        flux_err = get_value(flux_err, mask)
        
        # If the magnitude is not available but we have the flux, calculate the AB magnitude
        # TODO: add a unit for flux, now assuming nJy
        if mag is None and flux is not None:
            mag, mag_err = astro.nJy_to_ABmag(flux, flux_err)

        if mag is not None:
            # Apply extinction correction, if available
            if dered and ext is not None:
                mag = mag - ext
            elif dered and ext is None:
                logger.warning(f'Extinction correction is not available in catalog `{self.name}` for magnitude `{magnitude.filter}`.')

            return mag, mag_err
        
        # TODO: if there is a filter column given, only return those rows where
        #       the column matches the raw filter name and set everything else to NaN
        #       this requires adding a raw name to the magnitude class

        return None, None

    def get_magnitude_by_name(self, magnitude: Magnitude, observed=None, dered=None, magnitude_type=None, mask=None):
        """
        Return the values of a magnitude from the data frame.
        """

        # Ignore magnitude type for now

        observed = observed if observed is not None else True
        dered = dered if dered is not None else True
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

        # TODO: check if this function is ever used, otherwise remove it and keep
        #       calculate_flux_filters_bands only

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

    def calculate_flux_filters_bands(self, unit=u.nJy):
        for p, phot in self.photometry.items():
            for m, mag in phot.magnitudes.items():
                if mag.columns is not None:
                    self.__calculate_flux_filter(mag, unit=unit)
                else:
                    raise NotImplementedError()

    def __calculate_flux_filter(self, magnitude: Magnitude, unit=u.nJy):
        """
        Calculate flux in filters. This is necessary to fill in missing data for
        the design files.
        """

        # TODO: column names are inconsistent with the version of this function that
        #       works with photometry objects

        filter = magnitude.get_name()           # Internal name of the filter in the form of {phot}_{filt}
        columns = magnitude.columns             # Columns containing magnitudes and fluxes
        filter_value = columns['filter_value'] if 'filter_value' in columns else None
        band = columns['band'] if 'band' in columns else None

        if filter_value is not None:
            mask = self.__data[columns['filter']] == filter_value
        else:
            mask = None

        any_flux = None         # First available flux value for a filter
        any_flux_err = None     # First available flux error for a filter
        
        # Do two iteration so that the missing columns can be filled in from
        # any_flux found in earlier iterations
        for prefix in ['', 'psf_', 'fiber_', 'total_']:
            flux_key = prefix + 'flux'
            flux_err_key = prefix + 'flux_err'
            mag_key = prefix + 'mag'
            mag_err_key = prefix + 'mag_err'

            flux_col = columns[flux_key] if flux_key in columns else None
            flux_err_col = columns[flux_err_key] if flux_err_key in columns else None
            mag_col = columns[mag_key] if mag_key in columns else None
            mag_err_col = columns[mag_err_key] if mag_err_key in columns else None

            flux_col_canonical = f'{filter}_{prefix}flux'
            flux_err_col_canonical = f'{filter}_{prefix}flux_err'

            # Read flux column or calculate flux from the magnitude if necessary
            flux = None
            flux_err = None
            if flux_col is None:
                # Flux column is not defined, check if the magnitude is available
                # and calculate the flux from it
                logger.warning(f'Missing column `{flux_key}` for filter {filter} in target list `{self.name}`.')

                if mag_col is not None and not np.all(self.__data[mag_col].isna()):
                    # The magnitude is available to calculate the flux from
                    logger.info(f'Calculating `{flux_key}` from `{mag_key}` in filter `{filter}`.')

                    mag = self.__data[mag_col]
                    mag_err = self.__data[mag_err_col] if mag_err_col is not None else None

                    flux = (np.array(mag) * u.ABmag).to_value(unit)
                    if mag_err is not None:
                        flux_err = 0.4 * np.log(10) * flux * np.array(mag_err)
                else:
                    # No flux or magnitude is available, set to NaN
                    logger.warning(f'No flux or magnitude of type `{prefix}` available for filter `{filter}` in target list `{self.name}`.')
            else:
                # The flux column is available, use it directly                
                flux = self.__data[flux_col]
                if flux_err_col is not None:
                    flux_err = self.__data[flux_err_col]

            # Apply the mask
            if mask is not None:
                if flux is not None:
                    flux[~mask] = np.nan
                if flux_err is not None:
                    flux_err[~mask] = np.nan

            if flux is None:
                logger.warning(f'Flux column `{flux_key}` for filter `{filter}` in target list `{self.name}` is None.')
            elif np.all(np.isnan(flux)):
                logger.warning(f'Flux column `{flux_key}` for filter `{filter}` in target list `{self.name}` is all NaN.')

            # Save the flux to the target list with canonical column name
            if flux_col_canonical not in self.__data:
                self.__data[flux_col_canonical] = flux if flux is not None else np.nan
            if flux_err_col_canonical not in self.__data:
                self.__data[flux_err_col_canonical] = flux_err if flux_err is not None else np.nan

            # Register the new columns in the photometry description
            if magnitude.columns is None:
                magnitude.columns = {}
            if flux_key not in magnitude.columns:
                magnitude.columns[flux_key] = flux_col_canonical
            if flux_err_key not in magnitude.columns:
                magnitude.columns[flux_err_key] = flux_err_col_canonical

            if any_flux is None:
                any_flux = flux
                any_flux_err = flux_err
            else:
                m = np.isnan(any_flux)
                if any_flux is not None and flux is not None:
                    any_flux[m] = flux[m]
                if any_flux_err is not None and flux_err is not None:
                    any_flux_err[m] = flux_err[m]

        # Fill in any possible missing or Nan values from any_flux
        for prefix in ['', 'psf_', 'fiber_', 'total_']:
            flux_col_canonical = f'{filter}_{prefix}flux'
            flux_err_col_canonical = f'{filter}_{prefix}flux_err'

            m = self.__data[flux_col_canonical].isna()
            if any_flux is not None:
                self.__data.loc[m, flux_col_canonical] = any_flux[m]
            if any_flux_err is not None:
                self.__data.loc[m, flux_err_col_canonical] = any_flux_err[m]
                    
