from ..photometry import Photometry, Magnitude
from ..io import ObservationSerializer
from .instrument import Instrument

class PanSTARRS(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('ps', latex=r'\mathrm{PS}\,')
        p.append_magnitude(Magnitude(
            'g',
            latex='g',
            conversion=0.8e33,
            sky=1e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'r',
            latex='r',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'i',
            latex='i',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'z',
            latex='z',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        return p