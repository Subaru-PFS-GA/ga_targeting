from ..util import *
from .photometry import Photometry
from .magnitude import Magnitude

class Color():
    def __init__(self, magnitudes, orig=None):
        if not isinstance(orig, Color):
            self.__magnitudes = magnitudes
        else:
            self.__magnitudes = magnitudes or safe_deep_copy(orig.magnitudes)
            
        # TODO: what if we have a color from different photometric systems?
        self.__photometry = self.__magnitudes[0].photometry

        self._update()
        self._validate()

    def __repr__(self):
        return f'Color({self.__magnitudes[0]}, {self.__magnitudes[1]})'

    def _update(self):
        pass

    def _validate(self):
        if len(self.__magnitudes) != 2:
            raise ValueError('Exactly two magnitudes must be specified.')

        p = None
        for m in self.__magnitudes:
            if p is None:
                p = m.photometry
            elif m.photometry != p:
                    raise Exception('Photometric systems must match.')

    def __get_magnitudes(self) -> Magnitude:
        return ReadOnlyList(self.__magnitudes)

    magnitudes = property(__get_magnitudes)

    def __get_photometry(self) -> Photometry:
        return self.__photometry

    photometry = property(__get_photometry)

    def get_latex(self):
        return f'{self.__photometry.latex} {self.__magnitudes[0].latex} - {self.__magnitudes[1].latex}'