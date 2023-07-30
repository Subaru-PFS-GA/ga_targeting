import logging
import numpy as np

from ..data import Catalog
from ..photometry import Color
from .linearselection import LinearSelection

class ColorSelection(LinearSelection):
    def __init__(self, color, min, max, orig=None):
        super(ColorSelection, self).__init__([color], 1.0, min, max, orig=orig)