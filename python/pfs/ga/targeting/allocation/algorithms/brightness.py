import numpy as np

from ...data import Catalog
from ...photometry import *
from ..associations import Associations
from ..targetingalgorithm import TargetingAlgorithm

class Brightness(TargetingAlgorithm):
    def __init__(self, magnitude: Magnitude = None, reverse=None, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, Brightness):
            self.__magnitude = magnitude
            self.__reverse = reverse or False
        else:
            self.__magnitude = magnitude or orig.__magnitude
            self.__reverse = reverse or orig.__reverse

    def __get_magnitude(self):
        return self.__magnitude

    def __set_magnitude(self, value):
        self.__magnitude = value

    magnitude = property(__get_magnitude, __set_magnitude)

    def argsort_assoc(self, catalog: Catalog, assoc: Associations, observed=None, mask=None, reverse=None):
        """
        Given the fiber associations, pick the brightest object available for the fiber.
        """
        observed = observed if observed is not None else catalog.observed
        mask = mask if mask is not None else np.full(catalog.shape, True, dtype=bool)
        reverse = reverse or self.__reverse

        m, s_m = catalog.get_magnitude(self.__magnitude, observed=observed, mask=mask)
        m = np.array(m)

        # For each fiber, order targets by decreasing brightness. The first
        # object will be the most preferred target
        if reverse:
            s = np.s_[::-1]
        else:
            s = np.s_[()]

        sorting = [ np.argsort(m[a])[s] for a in assoc.assoc ]

        return sorting