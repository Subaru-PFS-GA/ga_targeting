import numpy as np
import pandas as pd
import h5py
from typing import Callable

from ..data import Catalog
from .observationwriter import ObservationWriter

class Hdf5ObservationWriter(ObservationWriter):
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, Hdf5ObservationWriter):
            pass
        else:
            pass

    def write(self, catalog: Catalog, filename: str, target, inst, mask=None, name=None, mode='w'):
        group_name = f'obs/{target}/{inst}'

        if isinstance(mask, Callable):
            df = mask(catalog.data)
        elif mask is not None:
            df = catalog.data[mask]
        else:
            df = catalog.data

        with h5py.File(filename, mode) as h:
            # Recreate group if exists
            if group_name in h:
                del h[group_name]
            g = h.create_group(group_name)
            for col in df.columns:
                g.create_dataset(col, data=np.array(df[col]))

