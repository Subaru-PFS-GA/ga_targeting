from .axis import Axis

class XYAxis(Axis):
    def __init__(self, coord, limits=None, orig=None):
        super(XYAxis, self).__init__(limits=limits, orig=orig)

        if not isinstance(orig, XYAxis):
            self.__coord = coord
        else:
            self.__coords = coord or orig.__coord

        self.label = self.__coord

        self._validate()

    def _validate(self):
        pass

    def __get_coord(self):
        return self.__coord

    coord = property(__get_coord)