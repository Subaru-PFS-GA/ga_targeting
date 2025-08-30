from pfs.ga.common.util import *

from ..data import Catalog
from ..probabilitymap import ProbabilityMap
from . import Selection

class ProbabilityCut(Selection):
    """
    Implements a cut based on probability using a binned map.
    """

    def __init__(self, pmap: ProbabilityMap, population_id, lp_member_limit: float, orig=None):
        super(ProbabilityCut, self).__init__(orig=orig)

        if not isinstance(orig, ProbabilityCut):
            self.__pmap = pmap
            self.__population_id = population_id
            self.__lp_member_limit = lp_member_limit
        else:
            self.__pmap = pmap or safe_deep_copy(orig.pmap)
            self.__population_id = population_id or orig.__population_id
            self.__lp_member_limit = lp_member_limit or orig.__lp_member_limit

        self._validate()

    def _validate(self):
        pass

    def apply(self, catalog:Catalog, population_id=None, lp_member_limit=None, observed=None, mask=None):
        population_id = population_id or self.__population_id
        lp_member_limit = lp_member_limit or self.__lp_member_limit
        observed = observed if observed is not None else catalog.observed
        
        lp_member, lp_member_mask = self.__pmap.lookup_lp_member(catalog, observed=observed)
        
        mask = mask if mask is not None else np.full_like(lp_member_mask, True)

        mask[lp_member_mask] &= (lp_member[..., population_id][lp_member_mask] > self.__lp_member_limit)
        mask[~lp_member_mask] = False

        return mask
        

