import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.cm import get_cmap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from sklearn.neighbors import KDTree
from ics.cobraOps.Bench import Bench

from ..util import *
from ..diagram import styles
from ..allocation.associations import Associations
from ..allocation.fiberallocator import FiberAllocator
from ..projection import Pointing
from .instrument import Instrument
from .subaruwfc import SubaruWFC

class SubaruPFI(Instrument, FiberAllocator):
    """
    Fiber allocations with the Subaru PFS PFI unit
    """

    # IDs of cobras in corners
    CORNERS = [[742, 796, 56, 2338, 2392, 1652, 1540, 1594, 854]]

    # IDs of cobras in cobra block corners
    BLOCKS = [[742, 796, 56, 2],
                [1598, 2338, 2392, 1652],
                [854, 800, 1540, 1594]]

    def __init__(self, projection=None, orig=None):
        Instrument.__init__(self, orig=orig)
        FiberAllocator.__init__(self, orig=orig)

        if not isinstance(orig, SubaruPFI):
            self.__projection = projection or SubaruWFC(Pointing(0, 0))
        else:
            self.__projection = projection or orig.__projection

        self.__bench = Bench(layout='full')

    def __get_bench(self):
        return self.__bench

    bench = property(__get_bench)

    def get_cobra_centers(self):
        centers = np.array([self.__bench.cobras.centers.real, self.__bench.cobras.centers.imag]).T
        radec, mask = self.__projection.pixel_to_world(centers)

        return centers, radec

    def find_associations(self, *coords, mask=None):
        ctype, coords = normalize_coords(*coords)
        mask = mask if mask is not None else ()

        xy, fov_mask = self.__projection.world_to_pixel(coords[mask])
    
        b = self.__bench
        centers = np.array([b.cobras.centers.real, b.cobras.centers.imag]).T

        if fov_mask.sum() > 10:
            # Use kD-tree to find targets within R_max and
            # filter out results that are within R_min
            tree = KDTree(xy[fov_mask], leaf_size=2)
            outer = tree.query_radius(centers, b.cobras.rMax)
            inner = tree.query_radius(centers, b.cobras.rMin)

            assoc = [ np.setdiff1d(i, o, assume_unique=True) for i, o in zip(outer, inner) ]

            # Assoc is now a list with size equal to the number of cobras.
            # Each item in an integer array with the list of matching indexes into
            # the coords[mask][fov_mask] array.

            # Update associations to index into coords[mask] instead of coords[mask][fov_mask]
            index = np.arange(coords[mask].shape[0], dtype=np.int64)[fov_mask]
            assoc = [ index[a] for a in assoc ]
        else:
            # Do a 1-1 matching instead?
            assoc = [ np.array([], dtype=np.int64) for i in range(len(b.cobras.centers))]
        
        return Associations(assoc)
    
    def __plot_corners(self, ax):
        pass

    def __get_outline(self, ids, res, native_frame=None):
        # Calculate corners
        xy = []
        for i in ids:
            xy.append((self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag))
        xy.append(xy[0])    # wrap up
        xy = np.array(xy)

        # Calculate arcs
        axy = []
        t = np.linspace(0, 1, res)
        for xy0, xy1 in zip(xy[:-1], xy[1:]):
            axy.append(xy0 + t[:, np.newaxis] * (xy1 - xy0)[np.newaxis, :])
        
        axy = np.concatenate(axy)
        if native_frame == 'world':
            coords, mask = self.__projection.pixel_to_world(axy)
        elif native_frame == 'pixel':
            coords = axy
            mask = np.full(axy.shape[:-1], True, dtype=np.bool)
        
        return coords, mask

    def plot_focal_plane(self, ax: plt.Axes, diagram, res=None, **kwargs):
        corners = kwargs.pop('corners', False)
        blocks = kwargs.pop('blocks', False)

        style = styles.solid_line(**kwargs)
        res = res if res is not None else 36
        native_frame = diagram._get_native_frame()

        def plot_outline(ids):
            for ii in ids:
                xy, mask = self.__get_outline(ii, res, native_frame=native_frame)
                diagram.plot(ax, xy, mask=mask, native_frame=native_frame, **style)

        if corners:
            plot_outline(self.CORNERS)
        
        if blocks:
            plot_outline(self.BLOCKS)           

    def plot_cobras(self, ax, diagram, data=None, cmap='viridis', vmin=None, vmax=None, **kwargs):
        """
        Plots the fibers color-coded by the data vector
        """

        cmap = get_cmap(cmap) if not isinstance(cmap, Colormap) else cmap
        native_frame = diagram._get_native_frame()

        if data is not None:
            style = styles.closed_circle(**kwargs)
            vmin, vmax = find_plot_limits(data, vmin=vmin, vmax=vmax)
            color = get_plot_normalized_color(cmap, data, vmin=vmin, vmax=vmax)
        else:
            style = styles.open_circle(**kwargs)
            vmin, vmax = None, None
            color = None

        for i in range(len(self.__bench.cobras.centers)):
            # TODO: we are cheating here and assume that the patrol region will appear as a cirle
            #       in the plots, regardless of the actual projection

            xy = (self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag)
            rr = (self.__bench.cobras.centers[i].real + self.__bench.cobras.rMax[i], self.__bench.cobras.centers[i].imag)
            xy = diagram.project_coords(ax, xy, native_frame='pixel')
            rr = diagram.project_coords(ax, rr, native_frame='pixel')
            r = rr[0] - xy[0]

            if isinstance(color, Iterable):
                style['facecolor'] = color[i]
            elif color is not None:
                style['facecolor'] = color

            c = Circle(xy, r, **style)
            ax.add_patch(c)

        return ScalarMappable(cmap=cmap, norm=Normalize(vmin, vmax))

    def plot_fiber_numbers(self, ax, **kwargs):

        # TODO: figure out what coordinates and projections to use

        for i in range(len(self.__bench.cobras.centers)):
            x, y = (self.__bench.cobras.centers[i].real, self.__bench.cobras.centers[i].imag)
            ax.text(x, y, str(i), ha='center', va='center', fontsize=5)

    def plot_corners(self, ax, **kwargs):
        pass