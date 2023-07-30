from ..util import *

from .booleanselection import BooleanSelection

class AndSelection(BooleanSelection):
    def _operation(self, mask1, mask2):
        return (mask1 & mask2)