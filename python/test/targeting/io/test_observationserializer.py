import os

from test_base import TestBase
from pfs.ga.targeting.instrument.subaruhsc import SubaruHSC
from pfs.ga.targeting.io import ObservationSerializer

class ObservationSerializerTest(TestBase):
    def test_read_hdf5(self):
        fn = '/datascope/subaru/data/cmdfit/catalog/Munoz+18/munoz.h5'
        r = ObservationSerializer()
        r.read(fn, 'obs/umi/cfht')

    def test_read_hsc(self):
        fn = '/datascope/subaru/data/cmdfit/dSph/ursaminor_tpall3e_g24.cat'
        r = SubaruHSC.text_observation_reader()
        obs = r.read(fn)

        pass