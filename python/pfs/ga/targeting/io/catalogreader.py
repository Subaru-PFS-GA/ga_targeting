from ..util import *
from ..photometry import Photometry

class CatalogReader():
    def __init__(self, orig=None):
        if not isinstance(orig, CatalogReader):
            self.__photometry = {}
            self.__column_mapping = {}
        else:
            self.__photometry = safe_deep_copy(orig.__photometry)
            self.__column_mapping = safe_deep_copy(orig.__column_mapping)

    def __get_photometry(self):
        return ReadOnlyDict(self.__photometry)

    photometry = property(__get_photometry)

    def append_photometry(self, photometry):
        self.__photometry[photometry.name] = photometry

    def __get_column_mapping(self):
        return self.__column_mapping

    def __set_column_mapping(self, value):
        self.__column_mapping = value

    column_mapping = property(__get_column_mapping, __set_column_mapping)

    def _create_catalog(self):
        raise NotImplementedError()

    def read(self, filename, **kwargs):
        raise NotImplementedError()