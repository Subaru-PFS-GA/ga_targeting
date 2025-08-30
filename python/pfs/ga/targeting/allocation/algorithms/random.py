import numpy as np

from pfs.ga.common.data import Catalog
from pfs.ga.common.photometry import *

from ..associations import Associations
from ..targetingalgorithm import TargetingAlgorithm

class Random(TargetingAlgorithm):
        def __init__(self, orig=None):
            super().__init__(orig=orig)

        def argsort_assoc(self, assoc: Associations, mask=None, reverse=None):
            """
            Given the fiber associations, sort the objects randomly.
            """

            ix = [ np.random.permutation(np.arange(a.size, dtype=int)) for a in assoc.assoc ]

            return ix