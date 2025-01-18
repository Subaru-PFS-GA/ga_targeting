import matplotlib
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from ..util import *
from . import styles

class Diagram():
    def __init__(self, axes, orig=None):
        if not isinstance(orig, Diagram):
            self.__axes = axes
            self.__datalim = Bbox.null()
        else:
            self.__axes = axes or safe_deep_copy(orig.__axes)
            self.__datalim = orig.__datalim

    def __get_axes(self):
        return ReadOnlyList(self.__axes)

    axes = property(__get_axes)

    def apply(self, ax: plt.Axes, **kwargs):

        # First apply the autoscale limits
        ax.dataLim = self.__datalim
        ax.autoscale()

        # Override limits with user-specified limits

        ax.set_xlim(self.__axes[0].limits)
        ax.set_xlabel(self.__axes[0].label)
        if self.__axes[0].invert and not ax.xaxis_inverted():
            ax.invert_xaxis()

        ax.set_ylim(self.__axes[1].limits)
        ax.set_ylabel(self.__axes[1].label)
        if self.__axes[1].invert and not ax.yaxis_inverted():
            ax.invert_yaxis()
        
    def scatter(self, ax: plt.Axes, x, y, mask=None, s=None, scalex=True, scaley=True, **kwargs):
        style = styles.tiny_dots_scatter(**kwargs)
        
        s = s if s is not None else np.s_[:]
        mask = mask if mask is not None else np.s_[:]
        
        # scatter distinguishes between color and c
        color = style.pop('color', None)
        if isinstance(color, np.ndarray):
            color = color[mask][s]

        c = style.pop('c', None)
        if isinstance(c, np.ndarray):
            c = c[mask][s]

        size = style.pop('size', style.pop('s', None))
        if isinstance(size, np.ndarray):
            size = size[mask][s]

        style = styles.sanitize_style(**style)

        l = ax.scatter(x[mask][s], y[mask][s], color=color, c=c, s=size, **style)

        bbox = l.get_datalim(ax.transData)
        [[x0, y0], [x1, y1]] = bbox.get_points()
        self.__update_datalim([x0, x1], [y0, y1], scalex, scaley)

        self.apply(ax)

        return l

    def plot(self, ax: plt.Axes, x, y, fmt=None, mask=None, s=None, scalex=True, scaley=True, **kwargs):
        style = styles.tiny_dots_plot(**kwargs)

        s = s if s is not None else np.s_[:]
        mask = mask if mask is not None else np.s_[:]

        args = (x[mask][s], y[mask][s])
        if fmt is not None:
            args += (fmt,)
        
        lines = ax.plot(*args, **styles.sanitize_style(**style))
        
        for l in lines:
            x, y = l.get_data()
            self.__update_datalim(x, y, scalex, scaley)

        self.apply(ax)

        return lines

    def imshow(self, ax: plt.Axes, img, extent=None, scalex=True, scaley=True, **kwargs):
        # TODO: default style?
        style = styles.sanitize_style(**kwargs)

        im = ax.imshow(img, extent=extent, **style)

        # TODO: update datalim
        # self.__update_datalim([extent[0], extent[1]], [extent[2], extent[3]], scalex, scaley)

        self.apply(ax)

        return im
    
    def add_patch(self, ax: plt.Axes, patch, scalex=True, scaley=True):
        p = ax.add_patch(patch)

        # TODO: update axis limits selectively
        #       check _update_patch_limits
        pass
    
    def __update_datalim(self, x, y, scalex, scaley):   
        if scalex:
            self.__datalim.update_from_data_x(x, ignore=False)

        if scaley:
            self.__datalim.update_from_data_y(y, ignore=False)