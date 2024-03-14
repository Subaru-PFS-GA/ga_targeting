from datetime import datetime
from astropy.time import Time

class Pointing():
    def __init__(self, ra, dec, posang, obs_time=None, exp_time=None, nvisits=1):
        """
        Create a new `Pointing` instance that describes the telescope's position
        the time of the observation and the number and length of exposures.
        """
        self.__ra = ra
        self.__dec = dec
        self.__posang = posang
        self.__obs_time = obs_time if obs_time is not None else Time(datetime.now())
        self.__exp_time = exp_time if exp_time is not None else 900
        self.__nvisits = nvisits

    def __get_ra(self):
        return self.__ra
    
    def __set_ra(self, value):
        self.__ra = value

    ra = property(__get_ra, __set_ra)
    
    def __get_dec(self):
        return self.__dec
    
    def __set_dec(self, value):
        self.__dec = value

    dec = property(__get_dec, __set_dec)

    def __get_posang(self):
        return self.__posang
    
    def __set_posang(self, value):
        self.__posang = value

    posang = property(__get_posang, __set_posang)

    def __get_obs_time(self):
        return self.__obs_time
    
    def __set_obs_time(self, value):
        self.__obs_time = value

    obs_time = property(__get_obs_time, __set_obs_time)

    def __get_exp_time(self):
        return self.__exp_time
    
    def __set_exp_time(self, value):
        self.__exp_time = value

    exp_time = property(__get_exp_time, __set_exp_time)

    def __get_nvisits(self):
        return self.__nvisits
    
    def __set_nvisits(self, value):
        self.__nvisits = value

    nvisits = property(__get_nvisits, __set_nvisits)