import os

from test_base import TestBase
from pfs.ga.targeting.io import PfsFluxStdReader

class PfsFluxStdReaderTest(TestBase):
    def test_read(self):
        fn = '/datascope/subaru/data/cmdfit/dSph/fluxstd_ursaminor.feather'
        r = PfsFluxStdReader()
        fluxstd = r.read(fn)

        pass