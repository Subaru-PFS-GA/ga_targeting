from .axis import Axis

class RaDecAxis(Axis):
    def __init__(self, coord, limits=None, invert=None, orig=None):
        super(RaDecAxis, self).__init__(limits=limits, invert=invert, orig=orig)

        if not isinstance(orig, RaDecAxis):
            self.__coord = coord
        else:
            self.__coords = coord or orig.__coord

        if self.__coord == 'RA':
            self.label = r'$\alpha$'
        elif self.__coord == 'Dec':
            self.label = r'$\delta$'
        else:
            raise NotImplementedError()

        self._validate()

    def _validate(self):
        pass

    def __get_coord(self):
        return self.__coord

    coord = property(__get_coord)