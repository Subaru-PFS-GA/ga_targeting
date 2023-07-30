import pandas as pd

from ..util import safe_deep_copy
from ..photometry import Photometry
from ..data import Observation
from .observationreader import ObservationReader

class TextObservationReader(ObservationReader):
    def __init__(self, orig=None):
        super(TextObservationReader, self).__init__(orig=orig)

        if not isinstance(orig, TextObservationReader):
            self.__filter = None
            self.__column_names = []
            self.__kwargs = {}
        else:
            self.__filter = orig.__filter
            self.__column_names = safe_deep_copy(orig.__column_names)
            self.__kwargs = safe_deep_copy(orig.__kwargs)

    def __get_column_names(self):
        return self.__column_names


    def __set_column_names(self, value):
        self.__column_names = value

    column_names = property(__get_column_names, __set_column_names)

    def __get_filter(self):
        return self.__filter

    def __set_filter(self, value):
        self.__filter = value

    filter = property(__get_filter, __set_filter)

    def __get_kwargs(self):
        return self.__kwargs

    def __set_kwargs(self, value):
        self.__kwargs = value

    kwargs = property(__get_kwargs, __set_kwargs)

    def read(self, filename, filter=None, names=None, **kwargs):
        if self.__kwargs is not None:
            kwargs.update(self.__kwargs)
        names = names or self.__column_names
        
        df = pd.read_csv(filename, names=names, **kwargs)

        # Rename columns according to the mapping
        df.rename(columns = self.column_mapping, inplace = True)

        if self.__filter is not None:
            df = df[self.__filter(df)]

        # Reindex
        df.reset_index(drop=True, inplace=True)

        obs = self._create_catalog()
        obs._set_data(df)

        return obs