from ..data import Observation
from .catalogreader import CatalogReader

class ObservationReader(CatalogReader):
    def __init__(self, orig=None):
        super(ObservationReader, self).__init__(orig=orig)

        if not isinstance(orig, ObservationReader):
            pass
        else:
            pass

    def _create_catalog(self, name=None):
        obs = Observation(name=name)
        for _, p in self.photometry.items():
            obs.append_photometry(p)
        return obs