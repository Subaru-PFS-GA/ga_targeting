import os
import numpy as np

from test_base import TestBase
from pfs.ga.targeting.data import HistogramND

class HistogramNDTest(TestBase):
    def test_find_extents(self):
        data = np.random.uniform(-1, 1, size=(100, 3))
        hist = HistogramND()
        hist._HistogramND__find_extents(data, mask=None, extents=None, bin_sizes=(0.1, 0.1, 0.1))

    def test_find_bins(self):
        data = np.random.uniform(-1, 1, size=(100, 3))
        hist = HistogramND()
        hist._HistogramND__find_extents(data, mask=None, extents=None, bin_sizes=(0.1, 0.1, 0.1))
        hist._HistogramND__find_bins(None, bin_sizes=(0.1, 0.1, 0.1))

    def test_create_2d(self):
        data = np.random.uniform(-1, 1, size=(100, 2))
        hist = HistogramND()
        hist.create(data)
        self.assertEqual((50, 50), hist.hist.shape)

    def test_create_3d(self):
        data = np.random.uniform(-1, 1, size=(100, 3, 2))
        hist = HistogramND()
        hist.create(data)
        self.assertEqual((3, 50, 50), hist.hist.shape)

