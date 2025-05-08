from json import JSONDecoder

from .color import Color
from .magnitude import Magnitude
from .photometry import Photometry

class PhotometryDecoder(JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        if isinstance(obj, dict):
            if '__type' in obj:
                t = obj.pop('__type')
                if t == 'Photometry':
                    return Photometry(**obj)
                elif t == 'Color':
                    return Color(**obj)
                elif t == 'Magnitude':
                    return Magnitude(**obj)
                else:
                    raise ValueError(f'Unknown type: {t}')
        
        return obj