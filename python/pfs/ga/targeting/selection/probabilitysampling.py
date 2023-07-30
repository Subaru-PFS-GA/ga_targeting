from ..util import *
from ..data import Catalog
from ..probabilitymap import ProbabilityMap
from .selection import Selection

class ProbabilitySampling(Selection):
    """
    Implements a random selection that follows a probability map.
    """

    def __init__(self, pmap: ProbabilityMap, population_id, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, ProbabilitySampling):
            self.__pmap = pmap
            self.__population_id = population_id
        else:
            self.__pmap = pmap or safe_deep_copy(orig.pmap)
            self.__population_id = population_id or orig.__population_id

        self._validate()

    def _validate(self):
        pass

    def apply(self, catalog:Catalog, observed=None, population_id=None, mask=None):
        population_id = population_id or self.__population_id
        observed = observed if observed is not None else catalog.observed
        
        lp_member, lp_member_mask = self.__pmap.lookup_lp_member(catalog, observed=observed, mask=mask)
        lp = np.log(np.random.uniform(0, 1, size=lp_member[..., population_id].shape))
 
        mask = mask if mask is not None else np.full_like(lp_member_mask, True)
       
        mask[lp_member_mask] &= (lp_member[..., population_id][lp_member_mask] > lp[lp_member_mask])
        mask[~lp_member_mask] = False

        return mask