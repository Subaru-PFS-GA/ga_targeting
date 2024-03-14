from datetime import datetime

from ..util.args import *

# TODO: move this elsewhere
class Pointing():
    def __init__(self, *pos, posang=None, obs_time=None, orig=None):
        if not isinstance(orig, Pointing):
            self.__pos = normalize_pos(*pos)
            self.__posang = normalize_angle(posang, u.degree) or normalize_angle(0)
            self.__obs_time = normalize_time(obs_time or datetime.utcnow())
        else:
            self.__pos = normalize_pos(*pos) if pos is not None else orig.__pos
            self.__posang = normalize_angle(posang, u.degree) or orig.__posang
            self.__obs_time = obs_time or orig.__obs_time

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
