from ..photometry import Photometry, Magnitude
from ..io import ObservationSerializer
from .instrument import Instrument

class SubaruHSC(Instrument):

    @staticmethod
    def photometry():
        p = Photometry('hsc', latex=r'\mathrm{HSC}\,')
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
            'nb515',
            latex=r'\mathrm{NB}515',
            conversion=1e33,
            sky=1e5,
            sky_sigma=0,
            zero=0.09))
        return p

    @staticmethod
    def text_observation_reader(mags=None, ext=None):
        if mags is None:
            mags = ['i', 'g', 'nb515']
        
        if ext is None:
            ext = ['g', 'i', 'nb515']

        reader = ObservationSerializer(format='.csv')
        reader.append_photometry(SubaruHSC.photometry())
        reader.column_map = {
            'ID': 'objid',
            'RA': 'RA',
            'Dec': 'Dec',
        }
        reader.column_names = ['ID', 'RA', 'Dec', 'X', 'Y']

        for m in mags:
            reader.column_map[f'{m[0]}psf'] = f'obs_hsc_{m}'
            reader.column_map[f'{m[0]}psferr'] = f'err_hsc_{m}'

        for m in ext:
            reader.column_map[f'a_{m[0]}'] = f'ext_hsc_{m}'

        # These have to be done is separate loops because column order matters!

        for m in mags:
            reader.column_names.append(f'{m[0]}psf')

        for m in mags:
            reader.column_names.append(f'{m[0]}psferr')

        for m in mags:
            reader.column_names.append(f'cl{m[0]}')

        for m in ext:
            reader.column_names.append(f'a_{m[0]}')

        def filter(df):
            ff = None
            for m in mags:
                if m != 'n':
                    f = df[f'cl{m[0]}'] < 0.1
                    if ff is None:
                        ff = f
                    else:
                        ff = ff | f
            return ff

        reader.filter = filter
        reader.kwargs = dict(delimiter=r'\s+')

        return reader