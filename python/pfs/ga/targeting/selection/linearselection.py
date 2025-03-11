import logging
import numpy as np
from shapely.geometry import Point, Polygon

from ..util import *
from .util import *
from ..data import Catalog
from .selection import Selection

class LinearSelection(Selection):
    """
    Implements a straight cut selection in magnitudes or colors
    in the form of A x - b_min >= 0 and A x - b_max <= 0.
    """

    def __init__(self, axes, A, b_min=None, b_max=None, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, LinearSelection):
            self.__axes = axes
            self.__A = normalize_array(A, allow_none=False)
            self.__b_min = normalize_array(b_min)
            self.__b_max = normalize_array(b_max)
        else:
            self.__axes = axes or orig.__axes
            self.__A = normalize_array(A) or orig.__A
            self.__b_min = normalize_array(b_min) or orig.__b_min
            self.__b_max = normalize_array(b_max) or orig.__b_max

        self._validate()

    def _validate(self):
        validate_photometry(self.__axes)

    def __get_axes(self):
        return ReadOnlyList(self.__axes)

    axes = property(__get_axes)

    def get_shape(self,):
        if len(self.__axes) != 2:
            raise ValueError("Operation valid for 2D selections only.")

        corners = (
            [ self.__axes[0].limits[i] for i in (0, 1, 1, 0, 0)],
            [ self.__axes[1].limits[j] for j in (0, 0, 1, 1, 0)]
        )

        shape = Polygon([ Point(x, y) for x, y in zip(corners[0], corners[1]) ])

        for k, b in enumerate([self.__b_min, self.__b_max]):
            if b is not None:
                x = 2 * [None]
                y = 2 * [None]
                for i in range(self.__A.shape[1]):
                    for j in range(2):
                        x[j] = self.__axes[0].limits[j]
                        y[j] = (self.__b_min[i] - self.__A[0, i] * x[j]) / self.__A[1, i]

                # Now we have the two points to connect


        return shape

    def apply(self, catalog: Catalog, observed=None, mask=None):
        observed = observed if observed is not None else catalog.observed

        x = catalog.get_diagram_values(self.__axes, observed=observed)
        x = np.stack([i[0] for i in x], axis=0)

        # Rotate any additional dimensions (e.g. populations) to the beginning
        x = x.transpose(tuple(range(2, x.ndim)) + (0, 1))

        lhs = np.matmul(self.__A, x)

        # Rotate additional dimensions back to the end to match mask shape
        x = x.transpose((-2, -1) + tuple(range(x.ndim - 2)))
        lhs = lhs.transpose((-1,) + tuple(range(lhs.ndim - 1)))

        mask = mask.copy() if mask is not None else np.full_like(x[0], True, dtype=bool)
        
        if self.__b_min is not None:
            mask &= (lhs - self.__b_min >= 0)

        if self.__b_max is not None:
            mask &= (lhs - self.__b_max <= 0)

        return mask

