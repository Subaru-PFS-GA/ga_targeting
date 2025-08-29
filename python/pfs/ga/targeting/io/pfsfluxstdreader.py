import pandas as pd

from ..data import Observation
from . import ObservationSerializer

class PfsFluxStdReader(ObservationSerializer):
    """
    Reads sky coordinates from feather files
    """
    
    def __init__(self,
                 columns=None,
                 column_map=None,
                 orig=None,
                 **kwargs):
        
        super().__init__(
            columns = columns,
            column_map = column_map if column_map is not None else {
                'fluxstd_id': 'objid',
                'obj_id': 'orig_objid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax': 'parallax',
                'parallax_error': 'err_parallax',
                'pmra': 'pmra',
                'pmra_error': 'err_pmra',
                'pmdec': 'pmdec',
                'pmdec_error': 'err_pmdec',
            },
            orig=orig,
            **kwargs)

        if not isinstance(orig, PfsFluxStdReader):
            pass
        else:
            pass
