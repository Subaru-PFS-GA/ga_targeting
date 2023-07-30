import numpy as np
from scipy.optimize import curve_fit

from .radialprofile import RadialProfile

class Sersic(RadialProfile):
    def __init__(self, transformation=None, center=None, R_max=None, bins=None, orig=None):
        super(Sersic, self).__init__(transformation=transformation, center=center, R_max=R_max, bins=bins, orig=orig)

        if not isinstance(orig, Sersic):
            pass
        else:
            pass

    def __get_S_center(self):
        return self.__S(0, *self.params)

    S_center = property(__get_S_center)

    def __get_S_background(self):
        (S_0, S_e, R_e, n) = self.params
        return S_0

    S_background = property(__get_S_background)

    @staticmethod
    def __S(R, S_0, S_e, R_e, n=4):
        b_n = 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n**2) + 131 / (1148175 * n**3) - 2194697 / (30690717750 * n**4)
        return S_0 + S_e * np.exp(-b_n * ((R / R_e)**(1 / n) - 1))

    @staticmethod
    def __log_S(R, S_0, S_e, R_e, n=4):
        b_n = 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n**2) + 131 / (1148175 * n**3) - 2194697 / (30690717750 * n**4)
        return np.log(S_0 + S_e * np.exp(-b_n * ((R / R_e)**(1 / n) - 1)))

    def sample(self, size, *params, R_max=None):
        raise NotImplementedError()

    def prob(self, R, *params):
        raise NotImplementedError()

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
            p0 = [np.exp(log_S[-1]), np.exp(log_S[-1]) / 2, R[-1] / 2, 1.0]
        
        p, pcov = curve_fit(self.__log_S, R, log_S, sigma=log_S_sigma, p0=p0)
        self.params = p
        self.pcov = pcov
        return self.params, self.pcov

    def fit_nomorph(self, R=None, log_S=None, log_S_sigma=None):
        # Fit a restricted subset of parameters only

        if R is None or log_S is None:
            R, log_S, log_S_sigma = self.get_log_S()

        [S_0, S_e, R_e, n] = self.params

        def log_S_nomorph(R, S_0, S_e):
            return self.__log_S(R, S_0, S_e, R_e, n)

        p0 = [S_0, S_e]

        p, pcov = curve_fit(log_S_nomorph, R, log_S, sigma=log_S_sigma, p0=p0)
        self.params[:2] = p
        self.pcov = pcov
        return self.params, self.pcov