from datetime import datetime

class Pointing():
    def __init__(self, ra, dec, posang=None, time=None, orig=None):
        if not isinstance(orig, Pointing):
            self.__ra = ra or 0
            self.__dec = dec or 0
            self.__posang = posang or 0
            self.__time = time or datetime.utcnow()
        else:
            self.__ra = ra or orig.__ra
            self.__dec = dec or orig.__dec
            self.__posang = posang or orig.__posang
            self.__time = time or orig.__time

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

    def __get_time(self):
        return self.__time
    
    def __set_time(self, value):
        self.__time = value

    time = property(__get_time, __set_time)