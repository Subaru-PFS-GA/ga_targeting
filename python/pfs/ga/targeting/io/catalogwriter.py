from ..util import *
from ..photometry import Photometry

class CatalogWriter():
    def __init__(self, orig=None):
        if not isinstance(orig, CatalogWriter):
            self.__column_mapping = {}
        else:
            self.__column_mapping = safe_deep_copy(orig.__column_mapping)

    def __get_column_mapping(self):
        return self.__column_mapping

    def __set_column_mapping(self, value):
        self.__column_mapping = value

    column_mapping = property(__get_column_mapping, __set_column_mapping)

    def writer(self, catalog, filename, **kwargs):
        raise NotImplementedError()