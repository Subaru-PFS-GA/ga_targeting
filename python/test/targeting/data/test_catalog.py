import os
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

from test_base import TestBase
from pfs.ga.targeting.data import Catalog
from pfs.ga.targeting.io import GaiaReader

class CatalogTest(TestBase):
    def test_get_skycoords(self):
        gaia = GaiaReader().cone_search((0, 0), 20)
        coords = gaia.get_skycoords()
        coords = gaia.get_skycoords(include_pm=True)

        # Modify the epoch for a subset of the coordinates
        gaia.data.loc[np.s_[:100], 'epoch'] = 2000.0
        coords = gaia.get_skycoords(include_pm=True)

    def test_transform_coords(self):
        # Convert coordinates only
        gaia = GaiaReader().cone_search((0, 0), 20)
        gaia.transform_coords(target_frame='fk5', target_equinox=2015.0, apply_motion=False)
        
        # Apply spatial motion
        gaia = GaiaReader().cone_search((0, 0), 20)
        gaia.transform_coords(apply_motion=True, target_epoch=2025.0)

        gaia = GaiaReader().cone_search((0, 0), 20)
        gaia.epoch = None
        gaia.transform_coords(apply_motion=True, target_epoch=2025.0)
