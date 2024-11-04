import pandas as pd

from ..data import Observation
from .catalogserializer import CatalogSerializer
from .dataframeserializer import DataFrameSerializer

class ObservationSerializer(CatalogSerializer, DataFrameSerializer):
    def __init__(self, 
                 filters=None,
                 bands=None,
                 orig=None,
                 **kwargs):
        
        DataFrameSerializer.__init__(self, orig=orig, **kwargs)
        CatalogSerializer.__init__(self, filters=filters, bands=bands, orig=orig)

        if not isinstance(orig, ObservationSerializer):
            pass
        else:
            pass

    def _create_catalog(self, name=None):
        obs = Observation(name=name)
        for _, p in self.photometry.items():
            obs.append_photometry(p)
        return obs
    
    def read(self, filename: str, dataset=None, format=None, filter=None, **kwargs) -> pd.DataFrame:
        df = super().read(filename, dataset=dataset, format=format, mask=filter, **kwargs)
        
        obs = self._create_catalog()
        obs._set_data(df)

        return obs
    
    def write(self, catalog, filename: str, dataset=None, filter=None, **kwargs):
        super().write(catalog.data, filename, dataset=dataset, mask=filter, **kwargs)