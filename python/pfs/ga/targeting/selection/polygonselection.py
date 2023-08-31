import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
from shapely.geometry import Point, Polygon

from ..util import *
from .util import *
from ..data import Catalog
from .selection import Selection
from ..diagram import styles

class PolygonSelection(Selection):
    """
    Implements a 2D selection using polygons based on the
    shapely library
    """

    def __init__(self, axes, points, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, PolygonSelection):
            self.__axes = axes
            self.__shape = PolygonSelection._create_shape(points)
        else:
            self.__axes = axes or orig.__axes
            self.__shape = PolygonSelection._create_shape(points) or orig.__shape

        self._validate()

    def _validate(self):
        validate_photometry(self.__axes)

    @staticmethod
    def _create_shape(points):
        if points is None:
            return None
        else:
            return Polygon(points)

    def __get_axes(self):
        return ReadOnlyList(self.__axes)

    axes = property(__get_axes)

    def __set_points(self, points):
        self.__shape = PolygonSelection._create_shape(points)

    points = property(None, __set_points)

    def __get_shape(self) -> Polygon:
        return self.__shape

    def __set_shape(self, value: Polygon):
        self.__shape = value

    shape = property(__get_shape, __set_shape)

    def apply(self, catalog: Catalog, observed=None, mask=None):
        observed = observed if observed is not None else catalog.observed
        (x, _), (y, _) = catalog.get_diagram_values(self.__axes, observed=observed, mask=mask)
        mm = [ self.__shape.contains(Point(x, y)) for x, y in zip(np.ravel(x), np.ravel(y)) ]
        return np.array(mm).reshape(x.shape)

    def plot(self, ax: plt.Axes, **kwargs):
        style = styles.dashed_line(**kwargs)

        polygon = matplotlib.patches.Polygon(self.__shape.exterior.xy, True)
        ax.add_patch(polygon)