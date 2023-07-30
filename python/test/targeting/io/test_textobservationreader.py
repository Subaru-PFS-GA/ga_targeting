import os
from pfs.ga.targeting.instrument.subaruhsc import SubaruHSC

from test_base import TestBase
from pfs.ga.targeting.io import TextObservationReader
from pfs.ga.targeting.instrument import SubaruHSC

class TextObservationReaderTest(TestBase):
    def test_read(self):
        fn = '/datascope/subaru/data/cmdfit/dSph/umi_tpall3e_g24.cat'
        #r = TextObservationReader()
        r = SubaruHSC.text_observation_reader()
        r.read(fn)