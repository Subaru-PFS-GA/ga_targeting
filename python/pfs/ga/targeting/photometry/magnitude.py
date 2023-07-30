from ..util import *

class Magnitude():
    def __init__(self, filter=None, orig=None, **kwargs):
        conversion = kwargs.pop('conversion', None)
        sky = kwargs.pop('sky', None)
        sky_sigma = kwargs.pop('sky_sigma', None)
        zero = kwargs.pop('zero', None)
        softening = kwargs.pop('softening', None)

        if not isinstance(orig, Magnitude):
            self.__photometry = None
            self.__filter: str = filter
            self.__conversion: float = conversion or 1e6
            self.__sky: float = sky or 0
            self.__sky_sigma: float = sky_sigma or 0
            self.__zero: float = zero or 0
            self.__softening: float = softening or 0
        else:
            self.__photometry = orig.__photometry
            self.__filter: str = filter or safe_deep_copy(orig.__filter)
            self.__conversion: float = conversion or orig.__conversion
            self.__sky: float = sky or orig.__sky
            self.__sky_sigma: float = sky_sigma or orig.__sky_sigma
            self.__zero: float = zero or orig.__zero
            self.__softening: float = softening or orig.__softening

        self._validate()

    def _validate(self):
        pass

    def __get_photometry(self):
        return self.__photometry

    def _set_photometry(self, value):
        self.__photometry = value

    photometry = property(__get_photometry)

    def __get_filter(self) -> str:
        """
        Gets the filter id string.
        """

        return self.__filter

    def __set_filter(self, value: str):
        """
        Sets the filter id string.
        """

        self.__filter = value

    filter = property(__get_filter, __set_filter)

    def __get_conversion(self) -> float:
        """
        Gets the flux to counts conversion factor.
        """

        return self.__conversion

    def __set_conversion(self, value: float):
        """
        Sets the flux to counts conversion factor.
        """

        self.__conversion = value

    conversion = property(__get_conversion, __set_conversion)

    def __get_sky(self) -> float:
        """
        Gets the typical sky counts.
        """

        return self.__sky

    def __set_sky(self, value: float):
        """
        Sets the typical sky counts.
        """

        self.__sky = value

    sky = property(__get_sky, __set_sky)

    def __get_sky_error(self) -> float:
        """
        Gets the typical error of sky counts.
        """

        return self.__sky_sigma

    def __set_sky_error(self, value: float):
        """
        Sets the typical error of sky counts.
        """

        self.__sky_sigma = value

    sky_error = property(__get_sky_error, __set_sky_error)

    def __get_zero(self) -> float:
        """
        Gets the magnitude zero point shift.
        """

        return self.__zero

    def __set_zero(self, value: float):
        """
        Sets the magnitude zero point shift.
        """

        self.__zero = value

    zero = property(__get_zero, __set_zero)

    def get_name(self, prefix=None, name_mappings=None):
        prefix = prefix or ''
        name = prefix + self.__photometry.name + '_' + self.__filter
        if name_mappings is not None and name in name_mappings:
            name = name_mappings[name]
        return name

    def mag_to_sigma(self, mag):
        """
        Calculate the typical photometric error sigma for a given magnitude.
        """
        return ABmag_to_sigma(mag, self.__conversion, self.__sky, self.__sky_sigma, self.__zero, softening=self.__softening)