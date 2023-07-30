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
    def text_observation_reader():
        reader = TextObservationReader()
        reader.append_photometry(SubaruHSC.photometry())
        reader.column_mapping = {
            'ID': 'objid',
            'RA': 'RA',
            'Dec': 'Dec',
            'ipsf': 'obs_hsc_i',
            'gpsf': 'obs_hsc_g',
            'npsf': 'obs_hsc_nb515',
            'ipsferr': 'err_hsc_i',
            'gpsferr': 'err_hsc_g',
            'npsferr': 'err_hsc_nb515',
            'a_g': 'ext_hsc_g',
            'a_i': 'ext_hsc_i',
            'a_n': 'ext_hsc_nb515',
        }

        reader.column_names = ['ID', 'RA', 'Dec', 'X', 'Y', 'ipsf', 'gpsf', 'npsf', 'ipsferr', 'gpsferr', 'npsferr',
                'cli', 'clg', 'cln', 'a_g', 'a_i', 'a_n']

        reader.filter = lambda df: df['cli'] < 0.1
        reader.kwargs = dict(delimiter=r'\s+', engine='python')

        return reader