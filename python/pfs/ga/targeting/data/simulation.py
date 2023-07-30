import numpy as np

from ..util import *
from .catalog import Catalog
from ..photometry import Photometry, Color, Magnitude

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

    def apply_categories(self, a, g=None):
        if g is None and 'g' in self.__data:
            g = self.__data['g']
            
        if g is not None:
            idx = g.reshape((g.shape[0],) + (len(a.shape) - 1) * (1,))
            return np.take_along_axis(a, idx, axis=-1)
        else:
            return a

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
