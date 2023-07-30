import pandas as pd
import h5py

from .observationreader import ObservationReader

class Hdf5ObservationReader(ObservationReader):
    def __init__(self, orig=None):
        super(Hdf5ObservationReader, self).__init__(orig=orig)

        if not isinstance(orig, Hdf5ObservationReader):
            pass
        else:
            pass

    def read(self, filename, target, inst, mask=None, name=None, **kwargs):
        with h5py.File(filename, 'r') as h:
            g = h['obs'][target][inst]
            d = { k: g[k][()] for k in g.keys() }
        df = pd.DataFrame(d)

        # Rename columns according to the mapping
        df.rename(columns = self.column_mapping, inplace = True)

        if mask:
            df = df[mask(df)]

        # Reindex
        df.reset_index(drop=True, inplace=True)

        obs = self._create_catalog(name=name)
        obs._set_data(df)

        return obs
