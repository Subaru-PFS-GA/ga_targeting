from collections import defaultdict
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

    def __get_cobra_count(self):
        return self.__bench.cobras.nCobras

    cobra_count = property(__get_cobra_count)

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
            mask = np.full(axy.shape[:-1], True, dtype=bool)
        
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

            # TODO: This projects cobras into FP coordinates so it won't work 
            #       when we draw a FOV plot with some projection defines on the axis

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

    #region Netflow support functions

    def nf_get_closest_dots(self):
        """
        For each cobra, return the list of focal plane black dot positions.
        
        Returns
        =======
        list : List of arrays of complex focal plane positions.
        """
        res = []
        for cidx in range(len(self.__bench.cobras.centers)):
            nb = self.__bench.getCobraNeighbors(cidx)
            res.append(self.__bench.blackDots.centers[[cidx] + list(nb)])

        return res
    
    def nf_get_visibility_and_elbow(self, fp_pos):
        """
        Calculate the visibility and the corresponding elbow position for each
        target of a single visit from the focal plane positions in ˙tpos˙.

        Parameters
        ==========

        tpos : numpy.ndarray
            Complex focal plane positions of each target.

        Returns
        =======
        dict : Dictionary of lists of (cobra ID, elbow position) pairs, keyed by target indices.
        """

        from ics.cobraOps.TargetGroup import TargetGroup
        from ics.cobraOps.TargetSelector import TargetSelector

        class DummyTargetSelector(TargetSelector):
            def run(self):
                return

            def selectTargets(self):
                return

        tgroup = TargetGroup(fp_pos)
        tselect = DummyTargetSelector(self.__bench, tgroup)
        tselect.calculateAccessibleTargets()
        tmp = tselect.accessibleTargetIndices   # shape: (cobras, targets), padded with -1
        elb = tselect.accessibleTargetElbows    # shape: (cobras, targets), padded with 0+0j
        res = defaultdict(list)

        for cbr in range(tmp.shape[0]):
            for i, tidx in enumerate(tmp[cbr, :]):
                if tidx >= 0:
                    res[tidx].append((cbr, elb[cbr, i]))

        return res
    
    def nf_get_colliding_pairs(self, fp_pos, vis_elbow, dist):
        """Return the list of target pairs that would cause fiber
        collision when observed by two neighboring fibers.
        
        Parameters
        ==========

        fp_pos : array
            Complex focal plane positions of the targets
        vis_elbow: dict
            Visibility, dict of (cobra ID, elbow position), keyed by target index.
        dist:
            Maximum distance causing a collision.

        Returns
        =======
        set : pairs of target indices that would cause fiber top collisions.
        """

        # Collect targets associated with each cobra, for each visit
        fp_pos = np.array(fp_pos)
        ivis = defaultdict(list)
        for tidx, thing in vis_elbow.items():
            for (cidx, _) in thing:
                ivis[cidx].append(tidx)

        pairs = set()
        for cidx, i1 in ivis.items():
            # Determine target indices visible by this cobra and its neighbors
            nb = self.__bench.getCobraNeighbors(cidx)
            i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
            i2 = np.concatenate((i1, i2))
            i2 = np.unique(i2).astype(int)
            d = np.abs(np.subtract.outer(fp_pos[i1], fp_pos[i2]))
            for m in range(d.shape[0]):
                for n in range(d.shape[1]):
                    if d[m][n] < dist:
                        if i1[m] < i2[n]:               # Only store pairs once
                            pairs.add((i1[m], i2[n]))
        return pairs
    
    def nf_get_colliding_elbows(self, fp_pos, vis_elbow, dist):
        """
        For each target-cobra pair, and the corresponding elbow position,
        return the list of other targets that are too close to the "upper arm" of
        the cobra.
        
        Parameters
        ==========

        fp_pos : array
            Complex focal plane positions of the targets
        vis_elbow: dict
            Visibility, dict of (cobra ID, elbow position), keyed by target index.
        dist:
            Maximum distance causing a collision.

        Returns
        =======
        dict : Dictionary of list of targets too close to the cobra indexed by
               all possible target-cobra pairs.
        """

        # TODO: speed up by vectorizing loops

        # vis contains the visible cobra indices and corresponding elbow positions
        # by target index. Invert this and build dictionaries indexed by cobra
        # indices that contain the lists of targets with corresponding elbow positions.
        ivis = defaultdict(list)
        epos = defaultdict(list)    # target_index, elbow position pairs keyed by cobra_index
        for tidx, cidx_elbow in vis_elbow.items():
            for (cidx, elbowpos) in cidx_elbow:
                ivis[cidx].append(tidx)
                epos[cidx].append((tidx, elbowpos))

        res = defaultdict(list)
        for cidx, tidx_elbow in epos.items():
            # Determine target indices visible by neighbors of this cobra
            nb = self.__bench.getCobraNeighbors(cidx)       # list of cobra_index
            tmp = [ epos[j] for j in nb if j in epos ]      # all targets visible by neighboring cobras
            if len(tmp) > 0:
                i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
                i2 = np.unique(i2).astype(int)

                # For each target visible by this cobra and the corresponding elbow
                # position, find all targets which are too close to the "upper arm"
                # of the cobra
                for tidx, elbowpos in tidx_elbow:
                    ebp = np.full(len(i2), elbowpos)
                    tp = np.full(len(i2), fp_pos[tidx])
                    ti2 = fp_pos[i2]
                    d = self.__bench.distancesToLineSegments(ti2, tp, ebp)
                    res[(cidx, tidx)] += list(i2[d < dist])

        return res

    #endregion