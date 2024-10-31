import astropy.units as u

from ..util.args import *
from ..projection import Pointing

class Target():
    def __init__(self,
                 ID, name,
                 pos,
                 DM=None, DM_err=None,
                 dist=None, dist_err=None,
                 parallax=None, parallax_err=None,
                 pm=None, pm_err=None,
                 RV=None, RV_err=None,
                 equinox='J2000',
                 frame='icrs',
                 scale='utc',
                 target_class=None):

        self.__ID = ID
        self.__name = name

        # Sky position, including pm and rv
        # SkyCoords doesn't support errors so store them separately
        self.__pos, self.__pm_err, self.__RV_err, _ = normalize_skycoord(pos,
                                                                         pm=pm, pm_err=pm_err,
                                                                         RV=RV, RV_err=RV_err,
                                                                         equinox=equinox,
                                                                         frame=frame,
                                                                         scale=scale)

        self.__distance, self.__dist_err, self.__DM_err, self.__parallax_err = normalize_distance(dist, dist_err,
                                                                                                  DM=DM, DM_err=DM_err,
                                                                                                  parallax=parallax, parallax_err=parallax_err)
        
        self.__target_class = target_class

    #region Properties

    def __get_ID(self):
        """ID of the object : str"""
        return self.__ID
    
    ID = property(__get_ID)

    def __get_name(self):
        """name of the object : str"""
        return self.__name
    
    name = property(__get_name)
        
    def __get_pos(self):
        return self.__pos
        
    pos = property(__get_pos)

    def __get_distance(self):
        return self.__distance
    
    distance = property(__get_distance)

    def __get_target_class(self):
        return self.__get_target_class
    
    target_class = property(__get_target_class)

    def __get_ra(self):
        return self.__pos.ra
    
    ra = property(__get_ra)

    def __get_dec(self):
        return self.__pos.dec
    
    dec = property(__get_dec)

    def __get_pm(self):
        return np.stack(self.pos.pm_ra_cosdec, self.pos.pm_dec)
    
    pm = property(__get_pm)

    def __get_pm_err(self):
        return self.__pm_err
    
    pm_err = property(__get_pm_err)

    def __get_pmra(self):
        return self.__pos.pm_ra_cosdec
    
    pmra = property(__get_pmra)
    
    def __get_pmdec(self):
        return self.__pos.pm_dec
    
    pmdec = property(__get_pmdec)

    def __get_pmra_err(self):
        return self.__pm_err[0]
    
    pmra_err = property(__get_pmra_err)

    def __get_pmdec_err(self):
        return self.__pm_err[1]
    
    pmdec_err = property(__get_pmdec_err)

    def __get_RV(self):
        return self.__pos.radial_velocity
    
    RV = property(__get_RV)

    def __get_RV_err(self):
        return self.__RV_err
    
    RV_err = property(__get_RV_err)

    def __get_dist(self):
        return self.__distance.value
    
    dist = property(__get_dist)

    def __get_dist_err(self):
        return self.__dist_err
    
    dist_err = property(__get_dist_err)

    def __get_DM(self):
        return self.__distance.distmod
    
    DM = property(__get_DM)

    def __get_DM_err(self):
        return self.__DM_err
    
    DM_err = property(__get_DM_err)

    def __get_parallax(self):
        return self.__distance.parallax
    
    parallax = property(__get_parallax)

    def __get_parallax_err(self):
        return self.__parallax_err
    
    parallax_err = property(__get_parallax_err)

    def __get_equinox(self):
        # TODO: Does it have to be a number or a string is fine?
        return self.__pos.equinox
    
    equinox = property(__get_equinox)

    def __get_target_class(self):
        return self.__target_class
    
    target_class = property(__get_target_class)

    def get_center_pointing(self, posang=None):
        return Pointing(self.__pos, posang=posang)