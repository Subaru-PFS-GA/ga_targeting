from ..setup_logger import logger
from ..util import *
from .projection import Projection

class Transformation():
    """
    When implemented in derived classes, provides functions to transform
    the coordinates. This is not a projection but a transformation of
    focal plane coordinates.
    """

    def __init__(self, projection=None, orig=None):
        if not isinstance(orig, Transformation):
            self.__projection = projection
        else:
            self.__projection = projection or orig.__projection

    def __get_projection(self) -> Projection:
        return self.__projection

    def __set_projection(self, value: Projection):
        self.__projection = value

    projection = property(__get_projection, __set_projection)

    def _world_to_pixel(self, *coords):
        ctype, coords = normalize_coords(*coords)
        if self.__projection is not None:
            xy, _ = self.__projection.world_to_pixel(coords)
            return ctype, xy
        else:
            logger.warning('No projection is set for pixel frame transformation.')
            return ctype, coords

    def _pixel_to_world(self, ctype, xy):
        if self.__projection is not None:
            radec, _ = self.__projection.pixel_to_world(xy)
            return denormalize_coords(ctype, radec)
        else:
            return denormalize_coords(ctype, xy)

    def apply(self, *coords):
        raise NotImplementedError()

    def reverse(self, *coords):
        raise NotImplementedError()