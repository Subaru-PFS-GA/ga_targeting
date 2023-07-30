import os

from test_base import TestBase
from pfs.ga.targeting.io import Hdf5ObservationReader

class Hdf5ObservationReaderTest(TestBase):
    def test_read(self):
        fn = '/datascope/subaru/data/cmdfit/catalog/Munoz+18/munoz.h5'
        r = Hdf5ObservationReader()
        r.read(fn, 'umi', 'cfht')