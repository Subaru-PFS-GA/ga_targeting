from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

from ..util.args import *
from ..data import Observation
from .catalogserializer import CatalogSerializer

class TwoMASSReader(CatalogSerializer):
    """
    Class to query the 2MASS catalog using astroquery.
    Retrieves J, H, and K magnitudes and their errors for a given sky position.
    """

    def __init__(self,
                 catalog_name='2mass',
                 limits=None,
                 mask=None,
                 orig=None,
                 **kwargs):

        # 2MASS filter definitions (central wavelengths in microns)
        filters = {
            "2mass_J": dict(
                mag = "2mass_J",
                mag_err = "err_2mass_J",
            ), 
            "2mass_H": dict(
                mag = "2mass_H",
                mag_err = "err_2mass_K",
            ),
            "2mass_K": dict(
                mag = "2mass_K",
                mag_err = "err_2mass_K",
            ),
        }

        CatalogSerializer.__init__(self,
                                   catalog_name=catalog_name,
                                   filters=filters,
                                   # bands=bands,
                                   # photometry=photometry,
                                   limits=limits,
                                   orig=orig)

    def _create_catalog(self):
        obs = Observation(name=self.catalog_name, frame='icrs', epoch='J2000.0')
        for _, p in self.photometry.items():
            obs.append_photometry(p)
        return obs

    def cone_search(self, pos, rad):
        """
        Query the 2MASS catalog around a given position.

        Parameters
        ----------
        pos: `SkyCoord`:
            Cone search center
        rad: `Angle` or `Number`:
            Cone search radius, in arc min if `Number`

        Returns
        -------
        astropy.table.Table or None
            Table of results with J, H, K magnitudes and their errors, or None if no results found.
        """

        pos = normalize_pos(pos)
        rad = normalize_angle(rad, u.arcmin)

        vizier = Vizier(
            columns=[
                "RAJ2000", "DEJ2000",
                "Jmag", "e_Jmag",
                "Hmag", "e_Hmag",
                "Kmag", "e_Kmag"
            ]
        )

        result = vizier.query_region(pos, radius=rad, catalog="II/246/out")
        if result is None or len(result) == 0:
            return None
        
        df = result[0].to_pandas()
        df.rename(
            inplace=True,
            columns={
                "RAJ2000": "ra",
                "DEJ2000": "dec",
                "Jmag": "2mass_J",
                "e_Jmag": "err_2mass_J",
                "Hmag": "2mass_H",
                "e_Hmag": "err_2mass_H",
                "Kmag": "2mass_K",
                "e_Kmag": "err_2mass_K"
            }
        )

        # Create a `targetid` column by copying the index
        df['targetid'] = df.index

        self._read_photometry(df)

        obs = self._create_catalog()
        obs._set_data(df)
        obs._set_photometry(self.photometry)

        # If limits are defined, apply them
        # TODO: this could be done from the SQL query
        if self.limits is not None:
            obs = obs.apply_magnitude_limits(self.limits)

        return obs