from pfs.ga.common.photometry import Photometry, Magnitude

from .instrument import Instrument

class PS1(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('ps1', latex=r'\mathrm{PS1}\,')
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
        p.append_magnitude(Magnitude(
            'y',
            latex='y',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        return p