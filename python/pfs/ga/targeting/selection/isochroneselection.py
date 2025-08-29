from scipy.interpolate import interp1d

from pfs.ga.common.util import *
from pfs.ga.common.diagram import ColorAxis, MagnitudeAxis
from pfs.ga.common.selection import Selection

from ..data import Catalog

class IsochroneSelection(Selection):
    def __init__(self, isochrone, axes, selection_axis, selection_direction, DM=20, error_sigma=None, orig=None):
        super().__init__(orig=orig)
        
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
        (iso_x, iso_x_err), (iso_y, iso_y_err) = self.__isochrone.get_diagram_values(self.__axes, error_sigma)

        # Collect catalog magnitudes
        (cat_x, cat_x_err), (cat_y, cat_y_err) = catalog.get_diagram_values(self.__axes, observed=observed)

        # TODO: optionally, take error from catalog
        if error_sigma is not None:
            if iso_x_err is not None:
                iso_x = iso_x + error_sigma[0] * iso_x_err
            if iso_y_err is not None:
                iso_y = iso_y + error_sigma[1] * iso_y_err

        # Determine interpolation axis, which is the opposite as
        # the selection axis
        if selection_axis == 0:
            pass
        elif selection_axis == 1:
            # TODO: swap x and y values
            raise NotImplementedError()
        else:
            raise ValueError('Invalid interpolation axis.')
        
        # Interpolate
        if mode == 'extend':
            # Extend isochrone with a straight line
            fill_value = (iso_y[0], iso_y[-1])
        else:
            # Do not extend
            fill_value=np.nan

        ip = interp1d(iso_y, iso_x, bounds_error=False, fill_value=fill_value)

        mask = mask.copy() if mask is not None else np.full_like(cat_x, True, dtype=bool)
        if selection_direction == '+':
            mask &= (cat_x >= ip(cat_y))
        elif selection_direction == '-':
            mask &= (cat_x <= ip(cat_y))
        else:
            raise ValueError('Invalid selection_direction, must be + or -')

        return mask


