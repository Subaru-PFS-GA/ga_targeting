import os
import numpy as np
from datetime import datetime
from astropy import units as u

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing

class PointingTest(TestBase):
    def test_init(self):
        p = Pointing(1, 2)
        self.assertEqual(1, p.ra)
        self.assertEqual(2, p.dec)
        self.assertEqual(0, p.posang)
        self.assertIsNotNone(p.obs_time)

        p = Pointing(1, 2, posang=3)
        self.assertEqual(1, p.ra)
        self.assertEqual(2, p.dec)
        self.assertEqual(3, p.posang)
        self.assertIsNotNone(p.obs_time)

        p = Pointing(1, 2, obs_time=datetime.now())
        self.assertEqual(1, p.ra)
        self.assertEqual(2, p.dec)
        self.assertEqual(0, p.posang)
        self.assertIsNotNone(p.obs_time)

        p = Pointing(1, 2, posang=3, obs_time=datetime.now())
        self.assertEqual(1, p.ra)
        self.assertEqual(2, p.dec)
        self.assertEqual(3, p.posang)
        self.assertIsNotNone(p.obs_time)