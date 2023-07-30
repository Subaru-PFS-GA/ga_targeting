import matplotlib
import matplotlib.pyplot as plt

from ..util import *
from . import styles

class Diagram():
    def __init__(self, axes, orig=None):
        if not isinstance(orig, Diagram):
            self.__axes = axes
        else:
            self.__axes = axes or safe_deep_copy(orig.__axes)

    def __get_axes(self):
        return ReadOnlyList(self.__axes)

    axes = property(__get_axes)

    def apply(self, ax, **kwargs):
        ax.set_xlim(self.__axes[0].limits)
        ax.set_xlabel(self.__axes[0].label)
        if self.__axes[0].invert and not ax.xaxis_inverted():
            ax.invert_xaxis()

        ax.set_ylim(self.__axes[1].limits)
        ax.set_ylabel(self.__axes[1].label)
        if self.__axes[1].invert and not ax.yaxis_inverted():
            ax.invert_yaxis()
        
    def scatter(self, ax: plt.Axes, x, y, mask=None, s=None, **kwargs):
        style = styles.tiny_dots_scatter(**kwargs)
        
        s = s if s is not None else np.s_[:]
        mask = mask if mask is not None else np.s_[:]
        
        color = style.pop('color', style.pop('c', None))
        if isinstance(color, np.ndarray):
            color = color[mask][s]

        size = style.pop('size', style.pop('s', None))
        if isinstance(size, np.ndarray):
            size = size[mask][s]

        style = styles.sanitize_style(**style)

        l = ax.scatter(x[mask][s], y[mask][s], c=color, s=size, **style)
        self.apply(ax)

        return l

    def plot(self, ax: plt.Axes, x, y, fmt=None, mask=None, s=None, **kwargs):
        style = styles.tiny_dots_plot(**kwargs)

        s = s if s is not None else np.s_[:]
        mask = mask if mask is not None else np.s_[:]

        args = (x[mask][s], y[mask][s])
        if fmt is not None:
            args += (fmt,)
        
        l = ax.plot(*args, **styles.sanitize_style(**style))
        self.apply(ax)

        return l

    def imshow(self, ax: plt.Axes, img, extent=None, **kwargs):
        # TODO: default style?
        style = styles.sanitize_style(**kwargs)

        l = ax.imshow(img, extent=extent, **style)
        self.apply(ax)

        return l