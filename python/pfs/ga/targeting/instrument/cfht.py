from pfs.ga.common.photometry import Photometry, Magnitude

from .instrument import Instrument

class CFHT(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('cfht', latex=r'\mathrm{CFHT}\,')
        p.append_magnitude(Magnitude(
            'g',
            latex='g',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'r',
            latex='r',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        return p