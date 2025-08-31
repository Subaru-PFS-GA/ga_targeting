import numpy as np
from collections.abc import Iterable
from scipy.special import logsumexp
from scipy.ndimage import maximum_filter
import matplotlib.pyplot as plt

from pfs.ga.common.util import *
from pfs.ga.common.diagram import MagnitudeDiagram, SpatialDiagram, styles
from .data import HistogramND, Simulation, Catalog

import h5py

class ProbabilityMap(HistogramND):
    def __init__(self, axes, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, ProbabilityMap):
            self.__axes = axes
            self.__population_weights = None
        else:
            self.__axes = axes or safe_deep_copy(orig.__axes)
            self.__population_weights = safe_deep_copy(self.__population_weights)

    def __get_axes(self):
        return ReadOnlyList(self.__axes)

    axes = property(__get_axes)

    def __get_population_weights(self):
        return self.__population_weights
    
    population_weights = property(__get_population_weights)

    def from_simulation(self,
                        sim: Simulation,
                        population_weights=None,
                        merge_list=None,
                        observed=None,
                        use_p_stars=False,
                        mask=None,
                        extents=None,
                        bins=None,
                        bin_sizes=None,
                        digits=None):
        
        """
        Calculate the probability map for each population
        """
        
        observed = observed if observed is not None else sim.observed

        x = sim.get_diagram_values(self.__axes, observed=observed)
        x = np.stack([i[0] for i in x], axis=-1)

        # Use the probability value for each star. When stellar parameters are
        # sampled from a uniform distribution, use these for weighting the pmap
        if use_p_stars:
            weights = np.exp(sim.data['lp_stars'])
        else:
            weights = None
       
        # Generate the histograms
        super().create(x, weights=weights,
                       mask=mask, extents=extents, bins=bins, bin_sizes=bin_sizes, digits=digits)
        
        # Get weights from argument or simulation
        if population_weights is not None:
            self.__population_weights = population_weights
        else:
            self.__population_weights = sim.data['w']

        if merge_list is not None:
            self.merge_populations(merge_list, population_weights=None)

    def merge_populations(self, merge_list, population_weights=None):
        # Argument `merge_list` should be a list of list of populations to merge        

        if population_weights is not None:
            w = population_weights
        else:
            w = self.__population_weights

        # Number of new populations after merge
        nn = len(merge_list)

        # New population weights
        nw = np.zeros((nn,), dtype=w.dtype)

        # New histograms
        nhist = np.empty((nn,) + self.hist.shape[1:], dtype=self.hist.dtype)

        for i in range(nn):
            nw[i] = np.sum(w[merge_list[i]])
            nhist[i] = np.sum(w[merge_list[i]][:, np.newaxis, np.newaxis] * self.hist[merge_list[i]], axis=0)

        self.__population_weights = nw

        self._set_hist(nhist)

    def maximum_filter(self, filter_size=3):
        # Run a maximum filter over the map to fill in holes

        if not isinstance(filter_size, Iterable):
            filter_size = len(self.__axes) * [ filter_size ]

        filter_size = tuple(filter_size)
        nhist = maximum_filter(self.hist, size=(1,) + filter_size, mode='nearest', origin=0)

        # Renormalize each filtered histogram
        self._set_hist(nhist)

    def get_nonzero_mask(self, populations=None):
        # Generate a mask where number counts are zero in any of the
        # populations listed in `populations`.

        if populations is None:
            populations = range(self.hist.shape[0])

        mask = np.full(self.hist.shape[1:], True, dtype=bool)

        for p in populations:
            mask &= (self.hist[p] > 0)

        return mask

    def get_norm(self, mask=None):
        # Normalize by the sum of counts in each population, optionally using
        # a mask

        if mask is None:
            mask = True

        # shape: (populations, color1, color2) -> (populations)
        norm = np.sum(np.where(mask, self.hist, 0), axis=(-1, -2))
        return norm

    def get_lp_member(self, population_weights=None):
        """
        Returns the probability map for each population weighted by the population
        weights.
        """

        if population_weights is not None:
            w = population_weights
        else:
            w = self.__population_weights

        w /= np.sum(w)
        
        # Although simulations usually consits of the same number of simulated
        # stars for each population, we have to normalize probability for each population
        # since a mask is usually applied when the histograms are created.
        map_mask = self.get_nonzero_mask()
        norm = self.get_norm(mask=map_mask)
        nhist = w[:, np.newaxis, np.newaxis] * self.hist / norm[:, np.newaxis, np.newaxis]

        # Normalize the weighted sum of the histograms in each cell to calculate
        # memberhip probabilities        
        
        lp_member = np.log(nhist / np.sum(nhist, axis=0))
        mask_member = ~np.isnan(lp_member)
        lp_member[~mask_member] = -np.inf
        
        return lp_member, mask_member

    def lookup_lp_member(self, catalog: Catalog, population_weights=None, observed=None, mask=None):
        """
        Returns the color-diagram-based membership probability of stars using
        different populations weights for each.
        """

        observed = observed if observed is not None else catalog.observed
        map_mask = self.get_nonzero_mask()
        norm = self.get_norm(mask=map_mask)

        # Get colors / magnitudes
        x = catalog.get_diagram_values(self.__axes, observed=observed, mask=mask)
        x = np.array([i[0] for i in x])

        # Look up map indices of stars for each axis
        ixx = []
        mask_member = np.full_like(x[0], True, dtype=bool)
        for i, ax in enumerate(self.__axes):
            ix = np.digitize(x[i], self.edges[i], right=False)
            
            # Make sure the point is inside the coverage of the map
            mask_member &= (self.extents[i, 0] <= x[i]) & (x[i] < self.extents[i, 1]) & \
                    (0 <= ix) & (ix < self.edges[i].size - 1)
            
            ixx.append(ix)

        # TODO: the catalog could have a column for population weights

        if population_weights is not None:
            w = population_weights
        else:
            w = self.__population_weights

        # Based on the shape of population weights, we either have them for
        # each star (and population) or for each population
        if len(w.shape) == 2:
            # We have population weights for each star, leave it as is
            w = w[mask_member].T
        else:
            # We have weights for each population only
            w = w[..., np.newaxis]

        # TODO: w can be longer than x, then the mask needs to be applied
      
        # Contruct the index arrays, the first one indexes the populations
        sp = (slice(None),)
        ss = tuple([ ixx[i][mask_member] for i, _ in enumerate(self.__axes)])

        # Look up probabilities and apply population weights
        # We do not assume that the histograms are already normalized to have an integral of 1
        # p.shape = (populations, catalog[mask])
        p = (self.hist / norm[:, np.newaxis, np.newaxis])[sp + ss] * w

        # Normalize across populations to get membership probability
        p /= np.sum(p, axis=0)

        # Mask out stars with nan probability (zero weight in all populations)
        mask_member[mask_member][np.any(np.isnan(p), axis=0)] = False
        
        lp_member = np.full(mask_member.shape + (p.shape[0],), np.nan, dtype=float)
        lp_member[mask_member] = np.log(p.T)

        return lp_member, mask_member
    
    def create_random_mask(self, catalog: Catalog, lp_member=None,
                           population_weights=None, observed=None, mask=None):
        """
        Generate a random mask based on population membership probability
        """

        if lp_member is None:
            lp_member, mask_member = self.lookup_lp_member(catalog, population_weights=population_weights, observed=observed, mask=mask)

        p = np.exp(lp_member)
        p = np.where(np.isnan(p), 0, p)
        r = np.random.binomial(1, p) == 1

        return r

    def save(self, filename):
        with h5py.File(filename, 'w') as f:
            g = f.create_group('pmap')
            self.save_items(g)

    def save_items(self, g):
        super().save_items(g)
        g.create_dataset('population_weights', data=self.__population_weights)

    def load(self, filename):
        with h5py.File(filename, 'r') as f:
            g = f['pmap']
            self.load_items(g)

    def load_items(self, g):
        super().load_items(g)
        self.__population_weights = g['population_weights'][:]

    def plot(self, ax: plt.Axes, diagram, *args, **kwargs):
        if isinstance(diagram, MagnitudeDiagram):
            return self._plot_magnitude(ax, diagram, *args, **kwargs)
        else:
            raise NotImplementedError()

    def _plot_magnitude(self, ax: plt.Axes, cmd, population_id=0, **kwargs):
        style = styles.histogram_imshow(**kwargs)

        lp_member, _ = self.get_lp_member()
        l = cmd.imshow(ax, lp_member[population_id].T / np.log(10), extent=self.extents.flatten(), **style)
        cmd.apply(ax)

        return l