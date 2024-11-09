import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

import pfs.ga.targeting
from pfs.ga.targeting.instrument import SubaruHSC, SubaruPFI
from pfs.ga.targeting.config.instrumentoptionsconfig import InstrumentOptionsConfig
from pfs.ga.targeting.config.netflowoptionsconfig import NetflowOptionsConfig

class NetflowOptionsConfigTest(TestBase):
    def test_create_cobra_instrument_labels(self):
        instrument_options = InstrumentOptionsConfig.default()
        pfi = SubaruPFI(instrument_options=instrument_options)

        config = NetflowOptionsConfig()
        cobra_location_labels = pfi.generate_cobra_location_labels(ntheta=6)
        cobra_instrument_labels = pfi.generate_cobra_instrument_labels(ngroups=8)
        
        for s in range(4):
            lab = cobra_instrument_labels[pfi.fiber_map.cobraId[(pfi.fiber_map.spectrographId == s + 1) & (pfi.fiber_map.cobraId != 65535)] - 1]
            self.assertEqual(lab.min(), s * 8)
            self.assertEqual(lab.max(), (s + 1) * 8 - 1)