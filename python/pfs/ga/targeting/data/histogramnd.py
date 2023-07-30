import numpy as np
from collections.abc import Iterable

from ..util import *

class HistogramND():
    def __init__(self, orig=None):
        if not isinstance(orig, HistogramND):
            self.__extents = None
            self.__hist = None
            self.__bins = None
            self.__edges = None
        else:
            self.__extents = safe_deep_copy(orig.__extents)
            self.__hist = safe_deep_copy(orig.__hist)
            self.__bins = safe_deep_copy(orig.__bins)
            self.__edges = safe_deep_copy(orig.__edges)

    def __get_extents(self):
        return self.__extents

    extents = property(__get_extents)

    def __get_hist(self):
        return self.__hist

    def _set_hist(self, value):
        self.__hist = value

    hist = property(__get_hist)

    def __get_bins(self):
        return self.__bins

    bins = property(__get_bins)

    def __get_edges(self):
        return self.__edges

    edges = property(__get_edges)

    def __find_extents(self, data: np.ndarray, mask, extents, bin_sizes):
        if extents is None:
            extents = []
            
            for i in range(data[mask].shape[-1]):
                if bin_sizes is not None:
                    d = bin_sizes[i]
                    mn = np.ceil(np.min(data[mask][..., i]) / d) * d
                    mx = np.floor(np.max(data[mask][..., i]) / d) * d
                else:
                    mn = np.min(data[mask][..., i])
                    mx = np.max(data[mask][..., i])

                extents.append([mn, mx])

            extents = np.stack(extents, axis=0)
        else:
            extents = normalize_array(extents)

            
        self.__extents = extents

    def __find_bins(self, bins, bin_sizes):
        if bins is None:
            bins = []
            for i in range(self.__extents.shape[0]):
                mn, mx = self.__extents[i]

                if bin_sizes is not None:
                    b = int(np.ceil((mx - mn) / bin_sizes[i]))
                else:
                    b = 50

                bins.append(b)

        self.__bins = bins

    def create(self, data, weights=None, mask=None, extents=None, bins=None, bin_sizes=None, digits=None):
        mask = mask if mask is not None else np.full(data.shape[:-1], True, dtype=bool)
        weights = weights if weights is not None else np.full(data.shape[:-1], 1.0)

        bindefs = 0
        if bin_sizes is not None:
            bindefs += 1
            if not isinstance(bin_sizes, Iterable):
                bin_sizes = data[mask].shape[-1] * [ bin_sizes, ]

        if digits is not None:
            bindefs +=1
            if not isinstance(digits, Iterable):
                digits = data[mask].shape[-1] * [ digits, ]
            bin_sizes = [ 10 ** (-d) for d in digits ]

        if bins is not None:
            bindefs += 1

        if bindefs > 1:
            raise ValueError('Only one of `bins`, `bins_sizes` and `digits` can be specified.')

        self.__find_extents(data, mask, extents, bin_sizes)
        self.__find_bins(bins, bin_sizes)

        # Histogram dimensions should match last dimension of data array
        if data.shape[-1] == 2:
            # Histogram batch shape (i.e. how many populations)
            shp = tuple([slice(s) for s in data.shape[1:-1]])
            if len(shp) > 0:
                # Process additional dimensions
                idx = np.mgrid[shp]
                for i, v in np.ndenumerate(idx[0]):
                    sl = (slice(None), v)
                    hist, xedges, yedges = np.histogram2d(data[sl][mask[sl]][..., 0], data[sl][mask[sl]][..., 1], bins=self.__bins, range=self.__extents, weights=weights[sl][mask[sl]])

                    if self.__hist is None:
                        self.__hist = np.empty(data.shape[1:2] + hist.shape, dtype=hist.dtype)
                    self.__hist[v] = hist

                    if self.__edges is None:
                        self.__edges = [xedges, yedges]
            else:
                # Just 2 dimensions
                hist, xedges, yedges = np.histogram2d(data[mask][..., 0], data[mask][..., 1], bins=self.__bins, range=self.__extents, weights=weights)
                self.__hist = hist
                self.__edges = [xedges, yedges]
        else:
            raise NotImplementedError()
