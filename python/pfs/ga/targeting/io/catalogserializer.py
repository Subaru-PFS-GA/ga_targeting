from pandas import DataFrame

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
    
        if self.__filters is not None:
            # Filter names are listed explicitly in the `filters` attribute.
            filters = self.__filters
        elif self.__bands is not None:
            # The filter names are in a dataframe column. Get the unique values of this
            # column for each band and try to separate the filters into several photometric
            # systems.
            
            # Collect all unique filter names
            # TODO: collect column names from config for each band
            filters = {}
            for b, band in self.__bands.items():
                for f in df[band['filter']].unique():
                    filters[f] = band
        else:
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

        # Find the length of the longest key in prefixes
        pre_ml = max(map(lambda x: len(x), prefixes.keys()))
        post_ml = max(map(lambda x: len(x), postfixes.keys()))

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