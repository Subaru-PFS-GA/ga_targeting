from pandas import DataFrame

from ..setup_logger import logger

from ..util import *
from ..photometry import Photometry, Magnitude

class CatalogSerializer():

    def __init__(self, 
                 catalog_name=None,
                 filters=None,
                 bands=None,
                 photometry=None,
                 limits=None,
                 orig=None):
        
        """
        Implements function to serialize and deserialize catalogs.
        """
    
        if not isinstance(orig, CatalogSerializer):
            self.__catalog_name = catalog_name
            self.__filters = filters
            self.__bands = bands
            self.__photometry = photometry if photometry is not None else {}
            self.__limits = limits
        else:
            self.__catalog_name = catalog_name if catalog_name is not None else orig.__catalog_name
            self.__filters = filters if filters is not None else safe_deep_copy(orig.__filters)
            self.__bands = bands if bands is not None else safe_deep_copy(orig.__bands)
            self.__photometry = photometry if photometry is not None else safe_deep_copy(orig.__photometry)
            self.__limits = limits if limits is not None else safe_deep_copy(orig.__limits)

    #region Properties

    def _can_read(self, filename=None, format=None):
        # TODO: implement this method
        return True

    def _can_write(self, catalog=None, filename=None, format=None):
        # TODO: implement this method
        return True
    
    def __get_catalog_name(self):
        return self.__catalog_name
    
    def __set_catalog_name(self, catalog_name):
        self.__catalog_name = catalog_name

    catalog_name = property(__get_catalog_name, __set_catalog_name)
    
    def __get_filters(self):
        return self.__filters
    
    def __set_filters(self, filters):
        self.__filters = filters

    filters = property(__get_filters, __set_filters)

    def __get_bands(self):
        return self.__bands
    
    def __set_bands(self, bands):
        self.__bands = bands

    bands = property(__get_bands, __set_bands)

    def __get_photometry(self):
        return ReadOnlyDict(self.__photometry)

    photometry = property(__get_photometry)

    def __get_limits(self):
        return self.__limits
    
    def __set_limits(self, limits):
        self.__limits = limits

    limits = property(__get_limits, __set_limits)

    #endregion

    def append_photometry(self, photometry):
        self.__photometry[photometry.name] = photometry

    def _read_photometry(self, df: DataFrame, name=None, latex=None):
        """
        Generate the photometry objects based on the columns defined in the
        `filters` and `bands` columns of the data frame.
        """

        def is_filter_good(f, filter):
            # Check if the filter has a column for any type of photometry
            for prefix in [ '', 'psf_', 'fiber_', 'total_' ]:
                for type in [ 'flux', 'mag' ]:
                    k = f'{prefix}{type}'
                    if k in filter:
                        if filter[k] not in df.columns:
                            raise ValueError(f'Column {filter[k]} not found in dataframe.')

                        # If there's at least one non-NaN value in the column, consider the filter good
                        if (~df[filter[k]].isna()).any():
                            return True
                            
            return False

        def is_band_filter_good(f, band):
            # Only include the filter if there are number available in the column
            # Ignore columns with all NaN values
            for prefix in [ '', 'psf_', 'fiber_', 'total_' ]:
                for type in [ 'flux', 'mag' ]:
                    k = f'{prefix}{type}'
                    if k in band:
                        if band[k] not in df.columns:
                            raise ValueError(f'Column {band[k]} not found in dataframe.')

                        if (~df[band[k]].isna()).any():
                            filters[f] = band
                            return True
                        
            return False

        filters = {}

        # Start with "one column per filter".
        if self.__filters is not None:
            # Filter names are listed explicitly in the `filters` attribute.
            for f, filter in self.__filters.items():
                if is_filter_good(f, filter):
                    filters[f] = filter
        
        # Update with the "one column per band", filter names in a column.
        if self.__bands is not None:
            # The filter names are in a dataframe column. Get the unique values of this
            # column for each band and try to separate the filters into several photometric
            # systems.
            
            # Collect all unique filter names for each band
            for b, band in self.__bands.items():
                for f in df[band['filter']].unique():
                    if is_band_filter_good(f, band):
                        band = band.copy()
                        band['filter_value'] = f
                        band['band'] = b
                        filters[f] = band
                                    
        if len(filters) == 0:
            # No photometry information available
            return None

        try:
            photometry_names = self.__guess_photometry(filters, default_name=name, default_latex=latex)
            for p in photometry_names:
                photometry = Photometry(name=p, latex=p)
                for m, config in photometry_names[p].items():
                    magnitude = Magnitude(filter=m, latex=m)
                    magnitude._set_columns(config)
                    photometry.append_magnitude(magnitude)
                self.append_photometry(photometry)
        except ValueError:
            logger.error('Error guessing photometry names from filters. ')
            pass

        return self.__photometry
        
    def __guess_photometry(self, filters, default_name='unk', default_latex='unk'):
        """
        Guess the photometric systems and filter names from the provided filter names.
        """

        # Split filter names along underscores and look for common prefixes and
        # postfixes. The longer repeating prefix/postfix will be the name of photometric
        # system while the shorter ones are the filter names.
        prefixes = {}
        postfixes = {}
        single_names = {}
        for f, config in filters.items():
            parts = f.split('_')
            if len(parts) == 1:
                single_names[f] = config
            elif len(parts) == 2:
                prefix = parts[0]
                postfix = parts[-1]
        
                if prefix not in prefixes:
                    prefixes[prefix] = {}

                if postfix not in prefixes[prefix]:
                    prefixes[prefix][postfix] = config

                if postfix not in postfixes:
                    postfixes[postfix] = {}

                if prefix not in postfixes[postfix]:
                    postfixes[postfix][prefix] = config
            else:
                raise ValueError('Filter names can contain zero or one underscore.')
            
        # See whether postfix or prefix is longer, because that will be the name of
        # the photometric systems

        # Find the length of the shortest key in prefixes/postfixes
        # It's likely that the shortest is a filter name and the longest is the
        # name of the photometric system
        pre_ml = min(map(lambda x: len(x), prefixes.keys()))
        post_ml = min(map(lambda x: len(x), postfixes.keys()))

        if pre_ml > post_ml:
            photometry_names = prefixes
        else:
            photometry_names = postfixes

        # Add single names as an unknown photometric system
        if len(single_names) > 0:
            photometry_names['unknown'] = single_names

        return photometry_names

    def _create_catalog(self):
        raise NotImplementedError()