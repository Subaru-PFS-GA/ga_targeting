from pfs.ga.common.photometry import Photometry, Magnitude

from .instrument import Instrument
from ..io import ObservationSerializer

class CFHT(Instrument):

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

    @staticmethod
    def text_observation_reader(mags=None, ext=None, delimiter=r'\s+', skiprows=1):
        if mags is None:
            mags = ['g', 'r']
        
        if ext is None:
            ext = ['g', 'r']

        reader = ObservationSerializer(format='.csv')
        reader.append_photometry(CFHT.photometry())
        reader.column_map = {
            'ID': 'objid',
            'RA': 'RA',
            'Dec': 'Dec',
        }
        reader.column_names = ['ID', 'RA', 'Dec']

        # for m in mags:
        #     reader.column_map[f'{m[0]}psf'] = f'obs_cfht_{m}'
        #     reader.column_map[f'{m[0]}psferr'] = f'err_cfht_{m}'

        # for m in ext:
        #     reader.column_map[f'a_{m[0]}'] = f'ext_hsc_{m}'

        # These have to be done is separate loops because column order matters!

        for m in mags:
            reader.column_names.append(f'obs_sdss_{m[0]}')

        for m in mags:
            reader.column_names.append(f'err_sdss_{m[0]}')

        # for m in mags:
        #     reader.column_names.append(f'cl{m[0]}')

        for m in ext:
            reader.column_names.append(f'ext_sdss_{m[0]}')

        # def filter(df):
        #     ff = None
        #     for m in mags:
        #         if m != 'n':
        #             f = df[f'cl{m[0]}'] < 0.1
        #             if ff is None:
        #                 ff = f
        #             else:
        #                 ff = ff | f
        #     return ff

        # reader.filter = filter
        reader.kwargs = dict(
            delimiter=delimiter,
            skiprows=skiprows
        )

        return reader