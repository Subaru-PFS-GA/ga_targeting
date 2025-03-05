from collections.abc import Iterable
import numpy as np

from ..util import *
from .transformation import Transformation

class Whitening(Transformation):
    """
    Implements an iterative principal axes transformation to the input data
    to find the axes and position angle of an overdensity inside an irregular
    observation footprint.

    Deprojection is determined using a second moment method combined with a sigma clipping
    technique to acount for the jagged, non-elliptic edges of the observed field.
    The algorithm is very simple and does not support masks, varying magnitude limits,
    filling factors, etc.
    """

    def __init__(self, projection=None, orig=None):
        super().__init__(projection=projection, orig=orig)

        if not isinstance(orig, Whitening):
            self.__S = None
            self.__Vh = None
            self.__W = None
            self.__inv_W = None
            self.__M = None
            self.__iteration_callback = None
        else:
            self.__S = safe_deep_copy(orig.__S)
            self.__Vh = safe_deep_copy(orig.__Vh)
            self.__W = safe_deep_copy(orig.__W)
            self.__inv_W = safe_deep_copy(orig.__inv_W)
            self.__M = safe_deep_copy(orig.__M)
            self.__iteration_callback = safe_deep_copy(orig.__iteration_callback)  

    def __get_S(self):
        return self.__S

    S = property(__get_S)

    def __get_Vh(self):
        return self.__Vh

    Vh = property(__get_Vh)
    
    def __get_W(self):
        return self.__W

    W = property(__get_W)

    def __get_M(self):
        return self.__M

    M = property(__get_M)

    def __get_iteration_callback(self):
        return self.__iteration_callback

    def __set_iteration_callback(self, value):
        self.__iteration_callback = value

    iteration_callback = property(__get_iteration_callback, __set_iteration_callback)

    def create(self, *coords, iterations=10, s_cut=None, callback=None):
        ctype, xy = self._world_to_pixel(*coords)

        callback = callback if callback is not None else self.__iteration_callback
            
        # Sigma cuts as a function of iteration number
        if s_cut is None:
            s_cut = (iterations // 2) * [2,] + (iterations - iterations // 2) * [2.8,]

        if not isinstance(s_cut, Iterable):
            s_cut = iterations * [s_cut]

        # Start with all data points included
        mask = np.full((xy.shape[0],), True)

        for iter in range(iterations):
            M = np.mean(xy[mask, :], axis=0)    
            C = np.matmul((xy[mask] - M).T, (xy[mask] - M))
            U, S, Vh = np.linalg.svd(C / (np.sum(mask) - 1))
            R = np.sqrt(S)
            
            # Whitening matrix
            W = np.matmul(Vh.T, np.diag(1.0 / R)).T
            
            # Whiten and throw away everything outside s_cut sigma
            white = np.matmul(W, (xy - M).T).T

            if callback is not None:
                callback(xy, white, mask, M, W, U, S, Vh)

            # Since we divide by the square root of the eigenvalues of the
            # covariance matrix, sigma cuts are measured in unit radii
            mask = (white[:, 0]**2 + white[:, 1]**2 < s_cut[iter]**2)

        self.__M = M
        self.__S = S
        self.__Vh = Vh
        self.__W = W
        self.__inv_W = np.linalg.inv(self.__W)      # TODO: calculate from U, S, Vh

        return mask

    def apply(self, *coords):
        ctype, xy = self._world_to_pixel(*coords)
        white = np.matmul(self.__W, (xy - self.__M).T).T
        return denormalize_coords(ctype, white)

    def reverse(self, *white):
        ctype, white = normalize_coords(*white)
        xy = np.matmul(self.__inv_W, white.T).T + self.__M
        return self._pixel_to_world(ctype, xy)