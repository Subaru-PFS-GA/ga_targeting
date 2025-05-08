from json import JSONEncoder

from ..util import ReadOnlyDict, ReadOnlyList
from .color import Color
from .magnitude import Magnitude
from .photometry import Photometry

class PhotometryEncoder(JSONEncoder):
    def default(self, o):
        if isinstance(o, (tuple, list, ReadOnlyList)):
            # If the object is a list, we need to encode each item in the list
            return [ self.default(item) for item in o ]
        elif isinstance(o, (dict, ReadOnlyDict)):
            # If the object is a dict, we need to encode each key-value pair
            return { self.default(key): self.default(value) for key, value in o.items() }
        elif isinstance(o, (str, int, float, bool, type(None))):
            return o
        elif isinstance(o, Photometry):
            return dict(
                __type = 'Photometry',
                name = o.name,
                latex = o.latex,
                magnitudes = dict(o.magnitudes.items()),
                colors = list(o.colors),
            )
        elif isinstance(o, Color):
            return dict(
                __type = 'Color',
                magnitudes = list(o.magnitudes)
            )
        elif isinstance(o, Magnitude):
            return dict(
                __type = 'Magnitude',
                filter = o.filter,
                latex = o.latex,
                conversion = o.conversion,
                sky = o.sky,
                sky_sigma = o.sky_sigma,
                zero = o.zero,
                softening = o.softening,
                columns = o.columns,
            )
        
        return super().default(o)
        
