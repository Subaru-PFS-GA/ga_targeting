import numpy as np
from scipy.interpolate import interp1d

from ..util import *
from ..projection.transformation import Transformation

class RadialProfile():
    """
    Implements a base class for radial profile models including profile fitting and
    various other calculations.

    Profiles are fitted to a histogram of elliptical radii. The results of the
    deprojection algorithm must be provided in a form of a transformation, see e.g.
    `pfs.ga.targeting.projection.Whitening` for a method of second moments. Note,
    that the deprojection algorithm assumes non-spherical 2D coordinate, either using
    some map projection or focal plane projection when applicable.
    """

    def __init__(self, transformation=None, center=None, R_max=None, bins=None, orig=None):
        if not isinstance(orig, RadialProfile):
            self.__transformation = transformation
            self.__R_max = R_max
            self.__bins = bins or 10
            self.__hist = None
            self.__center = center
            self.__params = None
            self.__pcov = None
        else:
            self.__transformation = transformation or orig.__transformation
            self.__R_max = R_max or safe_deep_copy(orig.__R_max)
            self.__bins = bins or safe_deep_copy(orig.__bins)
            self.__hist = safe_deep_copy(orig.__hist)
            self.__center = center if center is not None else safe_deep_copy(orig.__center)
            self.__params = safe_deep_copy(orig.params)
            self.__pcov = safe_deep_copy(orig.__pcov)

    def __get_transformation(self) -> Transformation:
        return self.__transformation

    def __set_transformation(self, value: Transformation):
        self.__transformation = value

    transformation = property(__get_transformation)

    def __get_R_max(self):
        return self.__R_max

    def __set_R_max(self, value):
        self._R_max = value

    R_max = property(__get_R_max, __set_R_max)

    def __get_bins(self):
        return self.__bins

    def __set_bins(self, value):
        self.__bins = value

    bins = property(__get_bins, __set_bins)

    def __get_hist(self):
        return self.__hist

    def __set_hist(self, value):
        self.__hist = value

    hist = property(__get_hist, __set_hist)

    def __get_center(self):
        return self.__center

    def __set_center(self, value):
        self.__center = value

    center = property(__get_center, __set_center)

    def __get_params(self):
        return self.__params

    def __set_params(self, value):
        self.__params = value

    params = property(__get_params, __set_params)

    def __get_pcov(self):
        return self.__pcov

    def __set_pcov(self, value):
        self.__pcov = value

    pcov = property(__get_pcov, __set_pcov)

    S_center = property()
    S_background = property()

    def get_R_halflight(self):
        """
        Return the half-light radius for any profile. Evaluate the fitted
        profile at the histogram bin centers and interpolate.
        """

        S_0 = self.S_background
        S_c = self.S_center - S_0

        R = 0.5 * (self.__bins[1:] + self.__bins[:-1])
        A = np.pi * (self.__bins[1:] ** 2 - self.__bins[:-1] ** 2)
        S = self.__hist - self.S_background * A

        cs = np.cumsum(np.where(S > 0, S, 0))
        ip = interp1d(cs, R)

        return ip(cs.max() / 2)

    def R(self, *coords, center=None, mask=None):
        """
        Return radii from coordinates. When projection and transformation are
        properly set these should be the elliptical radii.
        """

        R, _ = self.world_to_elliptic(*coords, center=center, mask=mask)
        return R

    def phi(self, *coords, center=None, mask=None):
        """
        Return elliptical angles from coordinates.
        """

        _, phi = self.world_to_elliptic(*coords, center=center, mask=mask)
        return phi

    def world_to_elliptic(self, *coords, center=None, mask=None):
        """
        Return elliptical radii and angled from world coordinates.
        Projection and transformation must be properly set.
        """

        _, coords = normalize_coords(*coords)
        mask = mask if mask is not None else np.s_[()]

        if self.__transformation is not None:
            xy = self.__transformation.apply(coords[mask])
        else:
            xy = coords[mask]

        center = center if center is not None else self.__center
        if center is not None:
            xy = xy - center

        R = np.sqrt(np.sum((xy) ** 2, axis=1))
        phi = np.arctan2(xy[..., 1], xy[..., 0])

        return R, phi

    def elliptic_to_world(self, R, phi, center=None, ctype=None, mask=None):
        """
        Return coordinates from elliptical radii and angles. When projection and
        transformation are properly set these should be real world coordinates.
        """

        mask = mask if mask is not None else np.s_[()]
        center = center if center is not None else self.__center

        xy = (R[mask] * np.array([np.cos(phi[mask]), np.sin(phi[mask])])).T
        if center is not None:
            xy = center + xy

        if self.__transformation is not None:
            coords = self.__transformation.reverse(xy)

        if ctype is not None:
            return denormalize_coords(ctype, coords)
        else:
            return coords
       
    def histogram(self, *coords, mask=None, center=None):
        # Calculate a radial histogram from whitened coordinates
        
        R = self.R(*coords, mask=mask, center=center)        

        if self.__R_max is None:
            R_max = R.max()
        else:
            R_max = self.__R_max

        if isinstance(self.__bins, np.ndarray):
            bins = self.__bins
        else:
            bins = np.linspace(0, R_max, (self.__bins or 10) + 1)

        self.__hist, self.__bins = np.histogram(R[R <= R_max], bins=bins)

        return self.get_log_S()

    def get_log_S(self):
        R = 0.5 * (self.__bins[1:] + self.__bins[:-1])
        area = np.pi * (self.__bins[1:] ** 2 - self.__bins[:-1] ** 2)
        log_S = np.log(self.__hist / area)
        log_S_sigma = 1 / np.sqrt(self.__hist)

        return R, log_S, log_S_sigma

    def fit(self, R=None, log_S=None, log_S_sigma=None):
        raise NotImplementedError()

    def get_ellipse_xy(self, R, res=None, ctype='a'):
        """
        Returns a circle in transformed coordinates.
        """
        res = res if res is not None else 360

        phi = np.linspace(0, 2 * np.pi, res)
        xy = np.stack([R * np.cos(phi), R * np.sin(phi)], axis=-1)

        return denormalize_coords(ctype, xy)


    def get_ellipse(self, R, res=None, ctype='a'):
        """
        Returns the world coordinates of the ellipse that is a circle in deprojected
        coordinates.
        """
        xy = self.get_ellipse_xy(R, res=res, ctype=ctype)

        if self.__transformation is None:
            return xy
        else:
            return self.__transformation.reverse(xy)

    def plot_histogram(self, ax, **kwargs):
        R, log_S, log_S_sigma = self.get_log_S()

        l = ax.errorbar(R, log_S, 3 * log_S_sigma, **kwargs)
        
        return l

    def plot_profile(self, ax, R_min=None, R_max=None, **kwargs):
        R_min = R_min or 0
        R_max = R_max or self.__R_max

        R = np.linspace(R_min, R_max, 200)
        log_S = self.log_eval(R)
        
        l = ax.plot(R, log_S, **kwargs)

        return l