from datetime import datetime
from typing import Any
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import CartesianRepresentation
from astropy.coordinates.matrix_utilities import rotation_matrix
import astropy.units as u

from ..util.args import *

# TODO: move this elsewhere
class Pointing():
    def __init__(self, *pos, posang=None,
                 obs_time=None, exp_time=None,
                 priority=None, stage=None,
                 nvisits=None, orig=None):
        """
        Create a new `Pointing` instance that describes the telescope's position
        the time of the observation and the number and length of exposures.
        """

        if not isinstance(orig, Pointing):
            self.__pos = normalize_pos(*pos)
            self.__posang = normalize_angle(posang, u.degree) if posang is not None else normalize_angle(0)
            self.__obs_time = normalize_time(obs_time if obs_time is not None else datetime.utcnow())
            self.__exp_time = normalize_exp_time(exp_time if exp_time is not None else 900)
            self.__priority = priority
            self.__stage = stage
            self.__nvisits = nvisits
        else:
            self.__pos = normalize_pos(*pos) if (pos is not None and pos != ()) else orig.__pos
            self.__posang = normalize_angle(posang, u.degree) if posang is not None else orig.__posang
            self.__obs_time = normalize_time(obs_time) if obs_time is not None else orig.__obs_time
            self.__exp_time = normalize_exp_time(exp_time) if exp_time is not None else orig.__exp_time
            self.__priority = priority if priority is not None else orig.__priority
            self.__stage = stage if stage is not None else orig.__stage
            self.__nvisits = nvisits if nvisits is not None else orig.__nvisits

    def __repr__(self):
        arglist = f'{self.ra:0.4f}, {self.dec:0.4f}, posang={self.posang:0.4f}'
        if self.__exp_time is not None:
            arglist += f', exp_time={self.__exp_time}'
        if self.__obs_time is not None:
            arglist += f', obs_time={self.__obs_time}'
        if self.__stage is not None:
            arglist += f', stage={self.__stage}'
        if self.__priority is not None:
            arglist += f', priority={self.__priority}'
        return f'Pointing({arglist})'

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

    def __get_priority(self):
        return self.__priority
    
    def __set_priority(self, value):
        self.__priority = value

    priority = property(__get_priority, __set_priority)

    def __get_stage(self):
        return self.__stage
    
    def __set_stage(self, value):
        self.__stage = value

    stage = property(__get_stage, __set_stage)

    def __get_nvisits(self):
        return self.__nvisits
    
    def __set_nvisits(self, value):
        self.__nvisits = value

    nvisits = property(__get_nvisits, __set_nvisits)

    @staticmethod
    def __get_full_rotation(sep, pa, ra, dec):
        r1 = rotation_matrix(sep * u.deg, 'z')  # separation of the pointing center from the center of the object
        r2 = rotation_matrix(pa * u.deg, 'x')   # position angle of the ellipse
        r3 = rotation_matrix(dec * u.deg, 'y')    # declination
        r4 = rotation_matrix(-ra * u.deg, 'z')    # right ascension

        return r4 @ r3 @ r2 @ r1

    @staticmethod
    def __calculate_center(sep, pa, ra, dec):
        r = Pointing.__get_full_rotation(sep, pa, ra, dec)
        c = SkyCoord(CartesianRepresentation(1, 0, 0).transform(r))
        return c.ra.value, c.dec.value
    
    @staticmethod
    def from_relative_pos(*pos, dir, sep, posang=None,
                          obs_time=None, exp_time=None,
                          priority=None, stage=None,
                          nvisits=None):
        """
        Calculate the pointing center relative to the `ra` and `dec` coordinates in the direction
        `dir` at a separation angle of `sep`.
        """

        pos = normalize_pos(*pos)
        pra, pdec = Pointing.__calculate_center(sep, dir, pos.ra.value, pos.dec.value)
        return Pointing(pra, pdec, posang=posang,
                        obs_time=obs_time, exp_time=exp_time,
                        priority=priority, stage=stage,
                        nvisits=nvisits)