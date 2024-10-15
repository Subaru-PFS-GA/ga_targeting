from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

from ..util.args import *
from ..data import Observation
from .catalogreader import CatalogReader

class GaiaReader(CatalogReader):
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, GaiaReader):
            pass
        else:
            pass

    def cone_search(self, pos, rad):
        """Run a cone search on the GAIA archive.
        
        Parameters:
        -----------
        
        :pos: `SkyCoord`: Cone search center
        :rad: `Angle` or `Number`: Cone search radius, in arc min if `Number`"""

        # TODO: bring out parameters like data release etc.

        pos = normalize_pos(pos)
        rad = normalize_angle(rad, u.arcmin)

        
        gbprp = ['G', 'BP', 'RP']
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1
        
        query = f"""
            SELECT
                source_id, ra, dec, parallax, parallax_error,
                pm, pmra, pmra_error, pmdec, pmdec_error,
                phot_g_mean_mag, phot_g_mean_flux_over_error,
                phot_bp_mean_mag, phot_bp_mean_flux_over_error,
                phot_rp_mean_mag, phot_rp_mean_flux_over_error
            FROM gaiadr3.gaia_source
            WHERE CONTAINS(POINT('ICRS', gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec),
                           CIRCLE('ICRS', {pos.ra.degree:0.8f}, {pos.dec.degree:0.8f}, {rad.degree:0.8f})) = 1;"""
        print(query)
        job = Gaia.launch_job_async(query, dump_to_file=False)
        gaia = job.get_results()
        df = gaia.to_pandas()
        df.rename(inplace=True,
            columns={
                'source_id': 'objid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax': 'parallax',
                'parallax_error': 'err_parallax',
                'pm': 'pm',
                'pmra': 'pmra',
                'pmra_error': 'err_pmra',
                'pmdec': 'pmdec',
                'pmdec_error': 'err_pmdec',
                'phot_g_mean_mag': 'obs_gaia_g',
                'phot_g_mean_flux_over_error': 'err_gaia_g',
                'phot_bp_mean_mag': 'obs_gaia_bp',
                'phot_bp_mean_flux_over_error': 'err_gaia_bp',
                'phot_rp_mean_mag': 'obs_gaia_rp',
                'phot_rp_mean_flux_over_error': 'err_gaia_rp'
            })

        c = Observation('gaia', frame='icrs', equinox='J2015')
        c._set_data(df)

        return c
