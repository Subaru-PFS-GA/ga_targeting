from ..util import *
from ..data import Catalog
from ..probabilitymap import ProbabilityMap
from .selection import Selection

class MaskSelection(Selection):
    """
    Select objects based on a boolean mask.
    """

    def __init__(self, mask, orig=None):
        super(MaskSelection, self).__init__(orig=orig)

        if not isinstance(orig, MaskSelection):
            self.__mask = mask
        else:
            self.__mask = mask if mask is not None else orig.__mask

        self._validate()

    def _validate(self):
        pass

    def apply(self, catalog:Catalog, population_id=None, lp_member_limit=None, observed=None, mask=None):
        if mask is not None:
            return mask & self.__mask
        else:
            return self.__mask
        

