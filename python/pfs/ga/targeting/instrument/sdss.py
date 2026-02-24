from pfs.ga.common.photometry import Photometry, Magnitude

from .instrument import Instrument

class SDSS(Instrument):

    # This is SDSS photometry but with CFHT Megacam errors!

    @staticmethod
    def photometry():
        p = Photometry('sdss', latex=r'\mathrm{SDSS}\,')
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