from .axis import Axis

class ColorAxis(Axis):
    def __init__(self, color=None, limits=None, orig=None):
        super().__init__(limits=limits, orig=orig)

        if not isinstance(orig, ColorAxis):
            self.__color = color
        else:
            self.__color = color or orig.__color

        self.label = f'${self.__color.get_latex()}$'

        self._validate()

    def __repr__(self):
        return f'ColorAxis({self.__color}, limits={self.limits})'

    def _validate(self):
        if self.__color is None:
            raise ValueError('The value of color must be set.')

    def __get_color(self):
        return self.__color

    color = property(__get_color)

    def __get_photometry(self):
        return self.__color.photometry

    photometry = property(__get_photometry)