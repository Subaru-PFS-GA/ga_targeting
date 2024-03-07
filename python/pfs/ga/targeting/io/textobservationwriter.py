import pandas as pd
from typing import Callable

from ..util import safe_deep_copy
from ..photometry import Photometry
from ..data import Observation
from .observationwriter import ObservationWriter

class TextObservationWriter(ObservationWriter):
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, TextObservationWriter):
            self.__kwargs = {}
        else:
            self.__kwargs = safe_deep_copy(orig.__kwargs)

    def __get_kwargs(self):
        return self.__kwargs

    def __set_kwargs(self, value):
        self.__kwargs = value

    kwargs = property(__get_kwargs, __set_kwargs)

    def read(self, catalog, filename, mask=None, **kwargs):
        if self.__kwargs is not None:
            kwargs.update(self.__kwargs)

        df: pd.DataFrame = None

        if isinstance(mask, Callable):
            df = mask(catalog.data)
        elif mask is not None:
            df = catalog.data[mask]
        else:
            df = catalog.data

        if self.column_mapping is not None:
            # Implement column name mapping
            raise NotImplementedError()

        df.to_csv(filename, **kwargs)
