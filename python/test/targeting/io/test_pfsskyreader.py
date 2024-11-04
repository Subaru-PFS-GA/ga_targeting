import os

from test_base import TestBase
from pfs.ga.targeting.io import PfsSkyReader

class PfsSkyReaderTest(TestBase):
    def test_read(self):
        fn = '/datascope/subaru/data/cmdfit/dSph/sky_ursaminor.feather'
        r = PfsSkyReader()
        sky = r.read(fn)