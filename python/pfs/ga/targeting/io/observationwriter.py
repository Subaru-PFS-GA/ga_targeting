from ..data import Observation
from .catalogwriter import CatalogWriter

class ObservationWriter(CatalogWriter):
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, ObservationWriter):
            pass
        else:
            pass
