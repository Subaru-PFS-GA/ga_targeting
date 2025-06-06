import os

from test_base import TestBase
from pfs.ga.targeting.io import TwoMASSReader

class TwoMASSReaderTest(TestBase):
    def test_cone_search(self):
        r = TwoMASSReader()
        obs = r.cone_search([10, 20], 2)

        self.assertEqual((6, 9), obs.data.shape)

        self.assertTrue('2mass' in obs.photometry)
        self.assertTrue('H' in obs.photometry['2mass'].magnitudes)
        self.assertTrue('J' in obs.photometry['2mass'].magnitudes)
        self.assertTrue('K' in obs.photometry['2mass'].magnitudes)