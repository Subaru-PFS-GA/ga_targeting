import numpy as np

class Axis():
    def __init__(self, limits=None, label=None, invert=None, orig=None):
        if not isinstance(orig, Axis):
            self.__limits = limits or (None, None)
            self.__label = label
            self.__invert = invert if invert is not None else False
        else:
            self.__limits = limits or orig.__limits
            self.__label = label or orig.__label
            self.__invert = invert if invert is not None else orig.invert

    def __get_limits(self):
        return self.__limits

    def __set_limits(self, value):
        self.__limits = value

    limits = property(__get_limits, __set_limits)

    def __get_label(self):
        return self.__label

    def __set_label(self, value):
        self.__label = value

    label = property(__get_label, __set_label)

    def __get_invert(self):
        return self.__invert

    def __set_invert(self, value):
        self.__invert = value

    invert = property(__get_invert, __set_invert)

    def apply_limits(self, data, mask=None):
        mask = mask if mask is not None else np.full_like(data, True, dtype=bool)
        if self.__limits is not None:
            if self.__limits[0] is not None:
                mask &= (self.__limits[0] <= data)
            if self.__limits[1] is not None:
                mask &= (data <= self.__limits[1])
        return mask