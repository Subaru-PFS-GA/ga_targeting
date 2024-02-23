import numpy as np

from ...data import Catalog
from ...photometry import *
from ..associations import Associations
from ..targetingalgorithm import TargetingAlgorithm

class Priority(TargetingAlgorithm):
    def __init__(self, reverse=None, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, Priority):
            self.__reverse = reverse or False
        else:
            self.__reverse = reverse or orig.__reverse

    def argsort_assoc(self, priority, assoc: Associations, mask=None, reverse=None):
        """
        Given the fiber associations, pick the object with highest priority.
        """
        mask = mask if mask is not None else np.full(priority.shape, True, dtype=bool)
        reverse = reverse or self.__reverse

        # For each fiber, order targets by decreasing probability. The first
        # object will be the most preferred target
        if reverse:
            s = np.s_[::-1]
        else:
            s = np.s_[()]
            

        ix = [ np.argsort(priority[mask][a])[s] for a in assoc.assoc ]

        return ix
