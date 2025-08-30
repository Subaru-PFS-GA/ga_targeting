import os
import numpy as np
import json as json

from pfs.ga.common.photometry import *

from test_base import TestBase

class PhotometryTest(TestBase):
    def test_init(self):
        p = Photometry.PS1()

    def test_json_encode(self):
        p = Photometry.PS1()
        s = json.dumps(p, indent=4, sort_keys=False, cls=PhotometryEncoder)

        p1 = json.loads(s, cls=PhotometryDecoder)
        
        self.assertEqual(p.name, p1.name)
        self.assertEqual(p.latex, p1.latex)
        self.assertEqual(p.magnitudes.keys(), p1.magnitudes.keys())
        self.assertEqual(len(p.colors), len(p1.colors))