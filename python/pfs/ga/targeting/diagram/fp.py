from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt

from ..util import normalize_coords, denormalize_coords
from . import styles
from .spatialdiagram import SpatialDiagram
from .xyaxis import XYAxis

class FP(SpatialDiagram):
    """
    A plot of the focal plane of some instrument assuming a projection, potentially beyond WCS.
    """

    def __init__(self, projection=None, orig=None):
        axes = [XYAxis('x'), XYAxis('y')]
        super().__init__(axes, projection=projection, orig=orig)

        self._validate()

    def _validate(self):
        pass

    def _get_native_frame(self, native_frame=None):
        return native_frame if native_frame is not None else 'pixel'
    
    def apply(self, ax):
        super().apply(ax)
        ax.set_aspect('equal', adjustable='datalim')