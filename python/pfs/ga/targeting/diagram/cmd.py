import matplotlib
import matplotlib.pyplot as plt

from ..util import *
from ..data import Observation, Simulation
from .diagram import Diagram
from . import styles

class CMD(Diagram):
    def __init__(self, axes, orig=None):
        super().__init__(axes, orig=orig)

        if not isinstance(orig, CMD):
            pass
        else:
            pass

        self._validate()

    def _validate(self):
        pass

    def plot_selection(self, ax: plt.Axes, selection, **kwargs):
        selection.plot(ax, **kwargs)

    def can_plot(self, catalog, observed=None):
        from ..isochrone import Isochrone

        if isinstance(catalog, Observation):
            observed = observed if observed is not None else catalog.observed
            return catalog.has_diagram_values(self.axes, observed=observed)
        elif isinstance(catalog, Simulation):
            observed = observed if observed is not None else True
            return catalog.has_diagram_values(self.axes, observed=observed)
        elif isinstance(catalog, Isochrone):
            return catalog.has_diagram_values(self.axes, observed=False)

    def plot_isochrone(self, ax: plt.Axes, isochrone, observed=False, error_sigma=None, **kwargs):
        style = styles.dashed_line(**kwargs)

        (x, x_err), (y, y_err) = isochrone.get_diagram_values(self.axes, observed=observed)
        
        if error_sigma is not None:
            if x_err is not None:
                x = x + error_sigma[0] * x_err
            if y_err is not None:
                y = y + error_sigma[1] * y_err

        return self.plot(ax, x, y, **style)

    def plot_catalog(self, ax: plt.Axes, catalog, observed=None, population_id=None, apply_categories=None, g=None, mask=None, s=None, **kwargs):
        if isinstance(catalog, Observation):
            return self.plot_observation(ax, catalog, observed=observed, mask=mask, s=s, **kwargs)
        elif isinstance(catalog, Simulation):
            observed = observed if observed is not None else True
            apply_categories = apply_categories if apply_categories is not None else True

            return self.plot_simulation(ax, catalog, population_id=population_id, apply_categories=apply_categories, g=g,
                observed=observed, mask=mask, s=s, **kwargs)
        else:
            raise NotImplementedError()

    def plot_simulation(self, ax: plt.Axes, catalog, population_id=None, apply_categories=False, g=None,
            observed=None, mask=None, s=None, **kwargs):

        if population_id is not None:
            p = np.s_[..., population_id]
        else:
            p = np.s_[:]

        observed = observed if observed is not None else catalog.observed
        mask = mask[p] if mask is not None else np.s_[:]

        (x, _), (y, _) = catalog.get_diagram_values(self.axes, observed=observed)
        if apply_categories:
            x = catalog.apply_categories(x, g=g)
            y = catalog.apply_categories(y, g=g)

        return self.scatter(ax, x[p], y[p], mask=mask, s=s, **kwargs)

    def plot_observation(self, ax: plt.Axes, catalog, 
            observed=None, mask=None, s=None, **kwargs):
        
        observed = observed if observed is not None else catalog.observed
        (x, _), (y, _) = catalog.get_diagram_values(self.axes, observed=observed, mask=mask)

        return self.scatter(ax, x, y, s=s, **kwargs)

    def plot_probability_map(self, ax: plt.Axes, pmap, population_id=0, **kwargs):
        style = styles.histogram_imshow(**kwargs)

        lp_member, _ = pmap.get_lp_member()
        l = self.imshow(ax, lp_member[population_id].T / np.log(10), extent=pmap.extents.flatten(), **style)
        self.apply(ax)

        return l