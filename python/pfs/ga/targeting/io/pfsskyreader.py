import pandas as pd

from ..data import Observation
from .observationserializer import ObservationSerializer

class PfsSkyReader(ObservationSerializer):
    """
    Reads sky coordinates from feather files
    """
    
    def __init__(self,
                 columns=None,
                 column_map=None,
                 orig=None,
                 **kwargs):

        super().__init__(
            columns = columns if columns is not None else ['sky_id', 'ra', 'dec', 'tract', 'patch'],
            column_map = column_map if column_map is not None else {
                'sky_id': 'skyid',
                'ra': 'RA',
                'dec': 'Dec'
            },
            orig=orig,
            **kwargs)

        if not isinstance(orig, PfsSkyReader):
            pass
        else:
            pass

