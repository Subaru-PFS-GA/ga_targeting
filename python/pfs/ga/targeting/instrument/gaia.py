from pfs.ga.common.photometry import Photometry, Magnitude

from .instrument import Instrument

class Gaia(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('gaia', latex=r'\mathrm{Gaia}\,')
        p.append_magnitude(Magnitude(
            'bp',
            latex='bp',
            conversion=0.8e33,
            sky=1e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'rp',
            latex='rp',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'g',
            latex='g',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        return p