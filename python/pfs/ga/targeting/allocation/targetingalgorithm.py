import numpy as np

from ..data import Catalog

class TargetingAlgorithm():
    def __init__(self, orig=None):
        if not isinstance(orig, TargetingAlgorithm):
            pass
        else:
            pass

    def argsort_assoc(self, catalog: Catalog, assoc, mask=None, reverse=None):
        raise NotImplementedError()

