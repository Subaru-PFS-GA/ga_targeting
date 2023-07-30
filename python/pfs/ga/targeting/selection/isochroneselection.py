from scipy.interpolate import interp1d

from ..util import *
from ..data import Catalog
from ..diagram import ColorAxis, MagnitudeAxis
from .selection import Selection

class IsochroneSelection(Selection):
    def __init__(self, isochrone, axes, selection_axis, selection_direction, DM=20, error_sigma=None, orig=None):
        super(IsochroneSelection, self).__init__(orig=orig)
        
        if not isinstance(orig, IsochroneSelection):
            self.__isochrone = isochrone
            self.__axes = axes
            self.__selection_axis = selection_axis
            self.__selection_direction = selection_direction
            self.__DM = DM
            self.__error_sigma = error_sigma or [1.0, 0.0]
        else:
            self.__isochrone = isochrone or safe_deep_copy(orig.__isochrone)
            self.__axes = axes or safe_deep_copy(orig.__axes)
            self.__selection_axis = selection_axis or orig.__selection_axis
            self.__selection_direction = selection_direction or orig.__selection_direction
            self.__DM = DM or orig.__DM
            self.__error_sigma = error_sigma or orig.__error_sigma

        self._validate()

    def _validate(self):
        if len(self.__axes) != 2:
            raise ValueError('Exactly two axes must be specified.')

    def apply(self, catalog: Catalog, observed=None, selection_axis=None, selection_direction=None, mode=None, error_sigma=None, mask=None):        
        observed = observed if observed is not None else catalog.observed
        selection_axis = selection_axis or self.__selection_axis
        selection_direction = selection_direction or self.__selection_direction
        error_sigma = error_sigma or self.__error_sigma
        
        # Calculate value pairs and add photometric error
        iso = self.__isochrone.get_blurred_values(self.__axes, error_sigma)

        # Collect catalog magnitudes
        x = catalog.get_diagram_values(self.__axes, observed=observed)

        # Determine interpolation axis, which is the opposite as
        # the selection axis
        if selection_axis == 0:
            a0, a1 = (1, 0)
        elif selection_axis == 1:
            a0, a1 = (0, 1)
        else:
            raise ValueError('Invalid interpolation axis.')

        # Interpolate
        if mode == 'extend':
            # Extend isochrone with a straight line
            ip = interp1d(iso[a0], iso[a1], bounds_error=False, fill_value=(iso[a0][0], iso[a0][-1]))
        else:
            # Do not extend
            ip = interp1d(iso[a0], iso[a1], bounds_error=False, fill_value=np.nan)

        mask = mask if mask is not None else np.full_like(x[a0][0], True, dtype=bool)
        if selection_direction == '+':
            mask &= (x[a1][0] >= ip(x[a0][0]))
        elif selection_direction == '-':
            mask &= (x[a1][0] <= ip(x[a0][0]))
        else:
            raise ValueError('Invalid selection_direction, must be + or -')

        return mask


