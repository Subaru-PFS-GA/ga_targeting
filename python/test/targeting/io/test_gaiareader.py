import os

from test_base import TestBase
from pfs.ga.targeting.io import GaiaReader

class Hdf5ObservationReaderTest(TestBase):
    def test_read(self):
        r = GaiaReader()
        obs = r.cone_search([10, 20], 1)

        self.assertEqual((4, 16), obs.data.shape)