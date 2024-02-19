from ..photometry import Photometry, Magnitude
from ..io import TextObservationReader
from .instrument import Instrument

class SubaruHSC(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('hsc')
        p.append_magnitude(Magnitude(
            'g',
            conversion=0.8e33,
            sky=1e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'r',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'i',
            conversion=0.8e33,
            sky=5e5,
            sky_sigma=0,
            zero=0))
        p.append_magnitude(Magnitude(
            'nb515',
            conversion=1e33,
            sky=1e5,
            sky_sigma=0,
            zero=0.09))
        return p

    @staticmethod
    def text_observation_reader(mags=None, ext=None):
        if mags is None:
            mags = ['g', 'i', 'n']
        
        if ext is None:
            ext = ['g', 'i', 'n']

        reader = TextObservationReader()
        reader.append_photometry(SubaruHSC.photometry())
        reader.column_mapping = {
            'ID': 'objid',
            'RA': 'RA',
            'Dec': 'Dec',
        }
        reader.column_names = ['ID', 'RA', 'Dec', 'X', 'Y']
                # , 'ipsf', 'gpsf', 'npsf', 'ipsferr', 'gpsferr', 'npsferr',
                #   'cli', 'clg', 'cln', 'a_g', 'a_i', 'a_n']

        for m in mags:
            reader.column_mapping[f'{m}psf'] = f'obs_hsc_{m}'
            reader.column_mapping[f'{m}psferr'] = f'err_hsc_{m}'

        for m in ext:
            reader.column_mapping[f'a_{m}'] = f'ext_hsc_{m}'

        # These have to be done is separate loops because column order matters!

        for m in mags:
            reader.column_names.append(f'{m}psf')

        for m in mags:
            reader.column_names.append(f'{m}psferr')

        for m in mags:
            reader.column_names.append(f'cl{m}')

        for m in ext:
            reader.column_names.append(f'a_{m}')

        def filter(df):
            ff = None
            for m in mags:
                if m != 'n':
                    f = df[f'cl{m}'] < 0.1
                    if ff is None:
                        ff = f
                    else:
                        ff = ff | f
            return ff

        reader.filter = filter
        reader.kwargs = dict(delimiter=r'\s+', engine='python')

        return reader