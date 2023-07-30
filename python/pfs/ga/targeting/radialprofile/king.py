import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy.coordinates import SkyCoord

from .radialprofile import RadialProfile

class King(RadialProfile):
    def __init__(self, transformation=None, center=None, R_max=None, bins=None, orig=None):
        super(King, self).__init__(transformation=transformation, center=center, R_max=R_max, bins=bins, orig=orig)

        if not isinstance(orig, King):
            pass
        else:
            pass

    def __get_S_center(self):
        return self.__S(0, *self.params)

    S_center = property(__get_S_center)

    def __get_S_background(self):
        (S_b, S_0, R_c, R_t) = self.params
        return S_b

    S_background = property(__get_S_background)

    @staticmethod
    def __S(R, S_b, S_0, R_c, R_t):
        A = (1 + R * R / (R_c * R_c)) ** (-0.5)
        B = (1 + R_t * R_t / (R_c * R_c)) ** (-0.5)
        return S_b + np.where(R < R_t, S_0 * (A - B) * (A - B), 0.0)

    @staticmethod
    def __log_S(R, S_b, S_0, R_c, R_t):
        return np.log(King.__S(R, S_b, S_0, R_c, R_t))

    @staticmethod
    def __S_norm(S_b, S_0, R_c, R_t):
        # Integral of S(R) from 0 to R_t
        A = R_c**2 + R_t**2
        B = R_t / R_c
        return S_b * R_t + S_0 * R_c * (R_c * R_t - 2 * R_c**2 * np.sqrt(1 + B**2) * np.arcsinh(B) + A * np.arctan(B)) / A

    @staticmethod
    def __RS(R, S_b, S_0, R_c, R_t):
        return R * King.__S(R, S_b, S_0, R_c, R_t)

    @staticmethod
    def __log_RS(R, S_b, S_0, R_c, R_t):
        return np.log(King.__RS(R, S_b, S_0, R_c, R_t))

    @staticmethod
    def __RS_norm(S_b, S_0, R_c, R_t):
        A = 1 + R_t**2 / R_c**2
        B = R_c**2 + R_t**2
        return 0.5 * (S_b * R_t**2 + R_c**2 * S_0 * (-3 + (R_c**2 * (-1 + 4 * np.sqrt(A))) / B + np.log(A)))

    def sample(self, size, *params, R_max=None):
        # The profile S(R) gives count / unit area so we should draw samples from R * S(R)
        # We do a simple rejection sampling with the constant function as majorant
        # The majorant is determined by estimating the maximum of R * S(R) numerically because
        # solving for it results in a complex expression

        if len(params) == 0:
            params = self.params

        (S_b, S_0, R_c, R_t) = params

        R_max = R_max or self.R_max or R_t

        R = np.linspace(0, R_max, 100)
        RS = R * self.__S(R, *params)
        RS_max = 1.1 * RS.max()

        R = np.empty(size, dtype=float)
        m = np.full_like(R, True, dtype=bool)
        
        # Rejection sampling
        ss = m.sum()
        while ss > 0:
            R[m] = np.random.uniform(0, R_max, size=ss)
            y = np.random.uniform(0, RS_max, size=ss)
            A = R[m] * self.__S(R[m], *params)
            m[m] = (y > A)
            ss = m.sum()

        return R

    def prob(self, R, *params):
        # Return S(R) normalized by its integral between 0 and R_t
        if len(params) == 0:
            params = self.params

        N = King.__S_norm(*params)
        S = King.__S(R, *params)
        return S / N
        
    def eval(self, R, *params):
        if len(params) == 0:
            params = self.params
        return self.__S(R, *params)

    def log_eval(self, R, *params):
        if len(params) == 0:
            params = self.params
        return self.__log_S(R, *params)

    def fit(self, R=None, log_S=None, log_S_sigma=None):
        # Fit all parameters including morphological

        if R is None or log_S is None:
            R, log_S, log_S_sigma = self.get_log_S()

        if self.params is not None:
            p0 = self.params
        else:
            p0 = [np.exp(log_S[-1]), np.exp(log_S[0]), R[-1] / 2, R[-1]]
        
        p, pcov = curve_fit(self.__log_S, R, log_S, sigma=log_S_sigma, p0=p0)
        self.params = p
        self.pcov = pcov
        return self.params, self.pcov

    def fit_nomorph(self, R=None, log_S=None, log_S_sigma=None):
        # Fit a restricted subset of parameters only

        if R is None or log_S is None:
            R, log_S, log_S_sigma = self.get_log_S()

        [S_b, S_0, R_c, R_t] = self.params

        def log_S_nomorph(R, S_b, S_0):
            return self.__log_S(R, S_b, S_0, R_c, R_t)

        p0 = [S_b, S_0]

        p, pcov = curve_fit(log_S_nomorph, R, log_S, sigma=log_S_sigma, p0=p0)
        self.params[:2] = p
        self.pcov = pcov
        return self.params, self.pcov

    def get_axes_world(self, R=1.0):
        """
        Return the length of the axes in radians.
        """

        R = np.array([0, R, R, R, R])
        phi = np.array([0, 0, np.pi / 2, np.pi, 3 / 2 * np.pi])

        ra, dec = self.elliptic_to_world(R, phi, ctype='t')

        c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
        l = c[0].separation(c[1:]).degree
        pa = c[0].position_angle(c[1:]).degree

        # Order axes by length, just in case, although the whitening transformation should take
        # care of it if the eigenvectors are properly ordered by decreasing eigenvalues.
        if l[0] > l[1]:
            a, b = l[0], l[1]
        else:
            a, b = l[1], l[0]
            pa = pa[[1, 0, 3, 2]]

        # Pick the smaller position angle
        if pa[0] < pa[2]:
            paa, pab = pa[0], pa[1]
        else:
            paa, pab = pa[2], pa[3]

        return a * 60.0, b * 60.0, paa      # arcmin, arcmin, degree

    def get_params_world(self):
        (S_b, S_0, R_c, R_t) = self.params
        R_t, _, _ = self.get_axes_world(R_t)
        R_c, _, _ = self.get_axes_world(R_c)

        return R_c, R_t