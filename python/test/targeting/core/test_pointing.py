import os
import numpy as np
import matplotlib.pyplot as plt

from ics.cobraOps.Bench import Bench

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing

class PointingTest(TestBase):
    def test_init(self):
        pointing = Pointing([0, 0], posang=15, exp_time=900)

        pointing = Pointing(0, 0, posang=15, exp_time=900)
        self.assertEqual(0.0, pointing.ra)
        self.assertEqual(0.0, pointing.dec)
        self.assertEqual(15.0, pointing.posang)
        self.assertEqual(900, pointing.exp_time)
