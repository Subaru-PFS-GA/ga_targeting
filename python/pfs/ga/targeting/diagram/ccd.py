from ..util import *
from .cmd import CMD

class CCD(CMD):
    def __init__(self, axes, orig=None):
        super(CCD, self).__init__(axes, orig=orig)

        if not isinstance(orig, CCD):
            pass
        else:
            pass

        self._validate()

    def _validate(self):
        pass