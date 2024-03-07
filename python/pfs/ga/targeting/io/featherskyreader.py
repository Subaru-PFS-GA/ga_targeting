import pandas as pd

from ..data import Observation
from .observationreader import ObservationReader

class FeatherSkyReader(ObservationReader):
    """
    Reads sky coordinates from feather files
    """
    
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, FeatherSkyReader):
            pass
        else:
            pass

    def read(self, filename):
        columns = ['sky_id', 'ra', 'dec', 'epoch', 'tract', 'patch']
        df = pd.read_feather(filename, columns)

        df.rename(inplace=True,
                  columns={
                      'sky_id': 'skyid',
                      'ra': 'RA',
                      'dec': 'Dec'
                  })
        
        c = self._create_catalog()
        c._set_data(df)

        return c
