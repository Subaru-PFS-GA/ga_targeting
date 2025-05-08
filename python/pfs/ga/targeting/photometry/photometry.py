from ..util import *
from .magnitude import Magnitude

class Photometry():
    def __init__(self, name=None, latex=None, magnitudes=None, colors=None, orig=None):
        if not isinstance(orig, Photometry):
            self.__name = name
            self.__latex = latex
            self.__magnitudes = magnitudes if magnitudes is not None else {}
            self.__colors = colors if colors is not None else []
        else:
            self.__name = name or orig.__name
            self.__latex = latex or orig.__latex
            self.__magnitudes = magnitudes if magnitudes is not None else safe_deep_copy(orig.__magnitudes)
            self.__colors = colors if colors is not None else safe_deep_copy(orig.__colors)

        self._update()
        self._validate()

    def _update(self):
        for magnitude in self.__magnitudes.values():
            magnitude._set_photometry(self)
        for color in self.__colors:
            for magnitude in color.magnitudes:
                magnitude._set_photometry(self)

    def _validate(self):
        pass

    # TODO: move to instruments?

    def SDSS():
        p = Photometry('sdss')
        p.append_magnitude(Magnitude('u'))
        p.append_magnitude(Magnitude('g'))
        p.append_magnitude(Magnitude('r'))
        p.append_magnitude(Magnitude('i'))
        p.append_magnitude(Magnitude('z'))
        return p

    def CFHT():
        p = Photometry('cfht')
        p.append_magnitude(Magnitude('u'))
        p.append_magnitude(Magnitude('g'))
        p.append_magnitude(Magnitude('r'))
        p.append_magnitude(Magnitude('i'))
        p.append_magnitude(Magnitude('z'))
        return p

    def PS1():
        p = Photometry('ps1')
        p.append_magnitude(Magnitude('g'))
        p.append_magnitude(Magnitude('r'))
        p.append_magnitude(Magnitude('i'))
        p.append_magnitude(Magnitude('z'))
        p.append_magnitude(Magnitude('y'))
        p.append_magnitude(Magnitude('w'))
        p.append_magnitude(Magnitude('n'))
        return p

    def copy(self):
        return Photometry(orig=self)

    def __get_name(self) -> str:
        return self.__name
    
    def __set_name(self, value: str):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_latex(self) -> str:
        """
        Gets filter name with latex formatting.
        """

        return self.__latex
    
    def __set_latex(self, value: str):
        """
        Sets the filter name with latex formatting.
        """

        self.__latex = value

    latex = property(__get_latex, __set_latex)

    def __get_magnitudes(self):
        return ReadOnlyDict(self.__magnitudes)

    magnitudes = property(__get_magnitudes)

    def append_magnitude(self, magnitude: Magnitude):
        self.__magnitudes[magnitude.filter] = magnitude
        magnitude._set_photometry(self)

    def __get_colors(self):
        return ReadOnlyList(self.__colors)

    colors = property(__get_colors)

    def append_color(self, color):
        self.__colors.append(color)
