import numpy as np
import matplotlib.pyplot as plt

from pfs.ga.common.data import Catalog
from pfs.ga.common.util import *
from pfs.ga.common.photometry import Photometry, Color, Magnitude

class Simulation(Catalog):
    def __init__(self, name=None, orig=None):
        super(Simulation, self).__init__(name=name, orig=orig)

        if not isinstance(orig, Simulation):
            self.__data = None
        else:
            self.__data = safe_deep_copy(orig.__data)

    def __len__(self):
        raise NotImplementedError()

    def __get_shape(self):
        return (self.__data['t'].shape)

    shape = property(__get_shape)

    def __get_data(self) -> np.ndarray:
        return self.__data

    def _set_data(self, data: np.ndarray):
        self.__data = data

    data = property(__get_data)

    def __get_observed(self):
        return False

    observed = property(__get_observed)

    def _get_observed_default(self, observed=None):
        return observed if observed is not None else True

    def apply_categories(self, a, g=None):
        if g is None and 'g' in self.__data:
            g = self.__data['g']
            
        if g is not None:
            idx = g.reshape((g.shape[0],) + (len(a.shape) - 1) * (1,))
            return np.take_along_axis(a, idx, axis=-1)
        else:
            return a
        
    def has_magnitude(self, magnitude: Magnitude, observed=False, dered=True):
        if observed:
            if magnitude.get_name('obs_') in self.__data:
                return True
        else:
            if magnitude.get_name() in self.__data:
                return True
                    
        return False

    def get_magnitude(self, magnitude: Magnitude, observed=False, dered=True, mask=None):
        mask = mask if mask is not None else slice(None)

        def get_value(k, mask):
            if k in self.__data:
                return self.__data[k][mask]
            else:
                return None

        # TODO: Currently no extinction in simulations

        mag = None
        err = None

        if observed:
            mag = get_value(magnitude.get_name('obs_'), mask)
        else: 
            mag = get_value(magnitude.get_name(), mask)
            
        if mag is None:
            raise Exception(f'Magnitude `{magnitude.get_name()}` not available for the simulation `{self.name}`')

        err = get_value(magnitude.get_name('err_'), mask)

        return mag, err

    def plot_cmd(self, ax: plt.Axes, diagram, population_id=None, apply_categories=False, g=None,
            observed=None, mask=None, s=None, **kwargs):

        observed = observed if observed is not None else self.observed
        mask = mask if mask is not None else np.s_[:]

        (x, _), (y, _) = self.get_diagram_values(diagram.axes, observed=observed)
        
        if apply_categories:
            x = self.apply_categories(x, g=g)
            y = self.apply_categories(y, g=g)

            if isinstance(mask, np.ndarray):
                mask = self.apply_categories(mask, g=g)

        if population_id is not None:
            p = np.s_[..., population_id]
            x = x[p]
            y = y[p]

            if isinstance(mask, np.ndarray):
                mask = mask[p]

        l =  diagram.scatter(ax, x, y, mask=mask, s=s, **kwargs)

        diagram.apply(ax)

        return l

    def plot_spatial(self, ax: plt.Axes, diagram, population_id=None, apply_categories=False, g=None,
            observed=None, mask_fov=False, mask=None, s=None, **kwargs):
        
        # TODO: implement
        raise NotImplementedError()