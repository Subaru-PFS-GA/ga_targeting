from datetime import datetime
from typing import Any
from astropy.time import Time

from ..util.args import *

# TODO: move this elsewhere
class Pointing():
    def __init__(self, *pos, posang=None, obs_time=None, exp_time=None, nvisits=None, orig=None):
        """
        Create a new `Pointing` instance that describes the telescope's position
        the time of the observation and the number and length of exposures.
        """

        if not isinstance(orig, Pointing):
            self.__pos = normalize_pos(*pos)
            self.__posang = normalize_angle(posang, u.degree) if posang is not None else normalize_angle(0)
            self.__obs_time = normalize_time(obs_time if obs_time is not None else datetime.utcnow())
            self.__exp_time = normalize_exp_time(exp_time if exp_time is not None else 900)
            self.__nvisits = nvisits
        else:
            self.__pos = normalize_pos(*pos) if (pos is not None and pos != ()) else orig.__pos
            self.__posang = normalize_angle(posang, u.degree) if posang is not None else orig.__posang
            self.__obs_time = normalize_time(obs_time) if obs_time is not None else orig.__obs_time
            self.__exp_time = normalize_exp_time(exp_time) if exp_time is not None else orig.__exp_time
            self.__nvisits = nvisits if nvisits is not None else orig.__nvisits

    def __get_pos(self):
        return self.__pos
    
    def __set_pos(self, value):
        self.__pos = normalize_pos(value)
    
    pos = property(__get_pos, __set_pos)

    def __get_ra(self):
        return self.__pos.ra.degree if self.__pos is not None else None
    
    ra = property(__get_ra)

    def __get_dec(self):
        return self.__pos.dec.degree if self.__pos is not None else None
    
    dec = property(__get_dec)

    def __get_posang(self):
        return self.__posang.degree
    
    def __set_posang(self, value):
        self.__posang = normalize_angle(value, u.degree)

    posang = property(__get_posang, __set_posang)

    def __get_obs_time(self):
        return self.__obs_time
    
    def __set_obs_time(self, value):
        self.__obs_time = normalize_time(value)

    obs_time = property(__get_obs_time, __set_obs_time)

    def __get_exp_time(self):
        return self.__exp_time
    
    def __set_exp_time(self, value):
        self.__exp_time = normalize_exp_time(value)

    exp_time = property(__get_exp_time, __set_exp_time)

    def __get_nvisits(self):
        return self.__nvisits
    
    def __set_nvisits(self, value):
        self.__nvisits = value

    nvisits = property(__get_nvisits, __set_nvisits)