from ..util import *

from .selection import Selection

class BooleanSelection(Selection):
    def __init__(self, selections=None, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, BooleanSelection):
            self.__selections = selections or []
        else:
            self.__selections = selections or orig.__selections

        self._validate()

    def _validate(self):
        pass

    def __get_selections(self):
        return self.__selections

    selections = property(__get_selections)

    def plot(self, cmd, ax, **kwargs):
        for s in self.__selections:
            s.plot(ax, cmd, **kwargs)

    def apply(self, catalog, observed=None, mask=None):
        observed = observed if observed is not None else catalog.observed
        
        for s in self.__selections:
            m = s.apply(catalog, observed=observed, mask=mask)
            mask = m if mask is None else self._operation(mask, m)
        
        return mask

    def _operation(self, mask1, mask2):
        raise NotImplementedError()