import numpy as np

from ..data import Catalog
from .pointing import Pointing

class Projection():
    def __init__(self, pointing=None, orig=None):
        if not isinstance(orig, Projection):
            self.__pointing = pointing or Pointing(0, 0)
        else:
            self.__pointing = pointing or orig.__pointing

    def __get_pointing(self):
        return self.__pointing

    def __set_pointing(self, value):
        self.__pointing = value

    pointing = property(__get_pointing, __set_pointing)

    def world_to_pixel(self, *coords, mask=None):
        raise NotImplementedError()

    def pixel_to_world(self, *coords, mask=None):
        raise NotImplementedError()
