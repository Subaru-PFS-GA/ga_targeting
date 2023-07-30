import logging
import numpy as np

from ..data import Catalog
from .linearselection import LinearSelection

class MagnitudeSelection(LinearSelection):
    def __init__(self, magnitude, min, max, orig=None):
        super(MagnitudeSelection, self).__init__([magnitude], 1.0, min, max, orig=orig)