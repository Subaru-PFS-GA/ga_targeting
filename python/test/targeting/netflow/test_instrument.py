import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

import pfs.ga.targeting
from pfs.ga.targeting.netflow import Instrument, Pointing

from test_base import TestBase

class InstrumentTest(TestBase):
    def test_init(self):
        inst = Instrument()
        inst = Instrument(instrument_options={})

    def test_get_grand_fiber_map(self):
        inst = Instrument()
        inst._Instrument__get_grand_fiber_map()

    def test_get_fiber_map(self):
        inst = Instrument()
        inst.get_fiber_map()

    def test_get_configured_bench(self):
        inst = Instrument()
        inst._Instrument__get_configured_bench()

    def test_get_bench(self):
        inst = Instrument()
        inst.get_bench()

        inst = Instrument(instrument_options={'layout': 'full'})
        inst.get_bench()

        inst = Instrument(instrument_options={'layout': 'calibration'})
        inst.get_bench()

    def test_radec_to_altaz(self):
        inst = Instrument()
        inst.radec_to_altaz(226.3, 67.5, 0.0, Time("2016-04-03T08:00:00Z"))

    def test_radec_to_fp_pos(self):
        inst = Instrument()
        pointing = Pointing(226.3, 67.5, 0, Time("2024-06-10T00:00:00.0Z"))
        inst.radec_to_fp_pos(pointing,
                             np.array(226.3),
                             np.array(67.5))