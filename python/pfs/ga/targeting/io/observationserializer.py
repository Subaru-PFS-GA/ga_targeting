import pandas as pd

from ..data import Observation
from .catalogserializer import CatalogSerializer
from .dataframeserializer import DataFrameSerializer
from ..photometry import Photometry, Magnitude

class ObservationSerializer(CatalogSerializer, DataFrameSerializer):
    def __init__(self,
                 catalog_name=None,
                 filters=None,
                 bands=None,
                 photometry=None,
                 limits=None,
                 mask=None,
                 orig=None,
                 **kwargs):
        
        DataFrameSerializer.__init__(self, orig=orig, mask=mask, **kwargs)
        CatalogSerializer.__init__(self,
                                   catalog_name=catalog_name,
                                   filters=filters,
                                   bands=bands,
                                   photometry=photometry,
                                   limits=limits,
                                   orig=orig)

        if not isinstance(orig, ObservationSerializer):
            pass
        else:
            pass

    def _create_catalog(self):
        obs = Observation(name=self.catalog_name)
        for _, p in self.photometry.items():
            obs.append_photometry(p)
        return obs
    
    def read(self, filename: str, dataset=None, format=None, filter=None, **kwargs) -> pd.DataFrame:       
        df = super().read(filename, dataset=dataset, format=format, mask=filter, **kwargs)
        self._read_photometry(df)
        
        obs = self._create_catalog()
        obs._set_data(df)
        obs._set_photometry(self.photometry)

        # If limits are defined, apply them
        if self.limits is not None:
            obs = obs.apply_magnitude_limits(self.limits)

        return obs
    
    def write(self, catalog, filename: str, dataset=None, filter=None, **kwargs):
        super().write(catalog.data, filename, dataset=dataset, mask=filter, **kwargs)