from ..util import *
from ..photometry import Photometry

class CatalogSerializer():

    def __init__(self, 
                 filters=None,
                 bands=None,
                 orig=None):
        
        """
        Implements function to serialize and deserialize catalogs.
        """
    
        if not isinstance(orig, CatalogSerializer):
            self.__filters = filters
            self.__bands = bands
            self.__photometry = {}
        else:
            self.__filters = safe_deep_copy(orig.__filters)
            self.__bands = safe_deep_copy(orig.__bands)
            self.__photometry = safe_deep_copy(orig.__photometry)

    #region Properties

    def _can_read(self, filename=None, format=None):
        # TODO: implement this method
        return True

    def _can_write(self, catalog=None, filename=None, format=None):
        # TODO: implement this method
        return True
    
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

    #endregion

    def append_photometry(self, photometry):
        self.__photometry[photometry.name] = photometry

    def _create_catalog(self):
        raise NotImplementedError()