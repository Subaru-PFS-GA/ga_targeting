from .axis import Axis

class MagnitudeAxis(Axis):
    def __init__(self, magnitude=None, limits=None, orig=None):
        super().__init__(limits=limits, orig=orig)

        if not isinstance(orig, MagnitudeAxis):
            self.__magnitude = magnitude
        else:
            self.__magnitude = magnitude or orig.__magnitude

        self.label = f'${self.__magnitude.get_latex()}$'
        self.invert = True

        self._validate()

    def __repr__(self):
        return f'MagnitudeAxis({self.__magnitude}, limits={self.limits})'

    def _validate(self):
        if self.__magnitude is None:
            raise ValueError('The value of magnitude must be set.')

    def __get_magnitude(self):
        return self.__magnitude

    magnitude = property(__get_magnitude)

    def __get_photometry(self):
        return self.__magnitude.photometry

    photometry = property(__get_photometry)