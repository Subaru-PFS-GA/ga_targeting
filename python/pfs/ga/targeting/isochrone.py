from collections.abc import Iterable

from pfs.ga.isochrones import tensorlib as tt
from pfs.ga.isochrones.isogrid import IsoGrid

from .util import *
from .photometry import Color, Magnitude
from .diagram.diagramvalueprovider import DiagramValueProvider

class Isochrone(DiagramValueProvider):
    def __init__(self, orig=None):
        if not isinstance(orig, Isochrone):
            self.__DM = None
            self.__magnitudes = None
            self.__values = None
            self.__eep = None
            self.__M_ini = None
        else:
            self.__DM = orig.__DM
            self.__magnitudes = safe_deep_copy(orig.__magnitudes)
            self.__values = safe_deep_copy(orig.__values)
            self.__eep = safe_deep_copy(orig.__eep)
            self.__M_ini = safe_deep_copy(orig.M_ini)

    def __get_values(self):
        return ReadOnlyDict(self.__values)

    values = property(__get_values)

    def __get_magnitudes(self):
        return ReadOnlyDict(self.__magnitudes)

    magnitudes = property(__get_magnitudes)

    def __get_M_ini(self):
        return self.__M_ini
    
    def __get_eep(self):
        return self.__eep
    
    EEP = property(__get_eep)
    
    M_ini = property(__get_M_ini)

    def from_isogrid(self, photometry, isogrid: IsoGrid, Fe_H, log_t, M_ini=None, EEP=None, DM=None, name_mappings=None, **kwargs):
        """
        Creates an isochrone instance from a an isochrone grid using the cmdfit library.
        """

        if M_ini is not None and EEP is not None:
            raise ValueError('Either M_ini or EEP can be defined.')

        if M_ini is None and EEP is None:
            M_ini = np.linspace(0.3, 1.2, 10000)

        if DM is not None:
            self.__DM = DM
        
        # TODO: we could do some fancy broadcasting here
        if M_ini is not None:
            iso_M_ini = tt.tensor(M_ini, dtype=tt.float64)
            iso_Fe_H = tt.full(iso_M_ini.shape, Fe_H, dtype=tt.float64)
            iso_log_t = tt.full(iso_M_ini.shape, log_t, dtype=tt.float64)
        elif EEP is not None:
            iso_EEP = tt.tensor(EEP, dtype=tt.float64)
            iso_Fe_H = tt.full(iso_EEP.shape, Fe_H, dtype=tt.float64)
            iso_log_t = tt.full(iso_EEP.shape, log_t, dtype=tt.float64)
        else:
            raise NotImplementedError()
        
        # Figure out what values we need to have all magnitudes of the photometric system
        # TODO: this could be made into a generic function
        name_mappings = name_mappings or {}
        values = []
        for i, m in enumerate(photometry.magnitudes.values()):
            name = m.get_name(name_mappings=name_mappings)
            if name in isogrid.values:
                values.append(isogrid.values[name])

        if EEP is not None:
            values.append(isogrid.M_ini)

        # Perform the isochrone interpolation
        if M_ini is not None:
            eep, _, iso_values, iso_mask = isogrid.interp3d(iso_Fe_H, iso_log_t, iso_M_ini, values=values)
        elif EEP is not None:
            iso_values = isogrid._interp3d_EEP(iso_Fe_H, iso_log_t, iso_EEP, values=values, update_index=True)
            eep = iso_EEP
            iso_mask = tt.full(eep.shape, False)
            M_ini = iso_values.pop(-1)
        else:
            raise NotImplementedError()
        
        self.__M_ini = M_ini[~iso_mask]
        self.__eep = np.array(eep)[~iso_mask]
        
        self.__values = {}
        self.__magnitudes = {}
        i = 0
        for m in photometry.magnitudes.values():
            name = m.get_name(name_mappings=name_mappings)
            if name in isogrid.values:
                self.__values[name] = np.array(iso_values[i][~iso_mask])
                self.__magnitudes[name] = m
                i += 1

    def has_magnitude(self, magnitude: Magnitude, observed=False, name_mappings=None):
        return magnitude.get_name(name_mappings=name_mappings) in self.__values
        
    def get_magnitude(self, magnitude: Magnitude, DM=None, observed=False, mask=None, name_mappings=None):
        DM = DM or self.__DM or 0
        m = self.__values[magnitude.get_name(name_mappings=name_mappings)] + DM
        s_m = magnitude.mag_to_sigma(m)

        return m, s_m

    def get_color(self, color: Color, observed=False, mask=None):
        m1, s_m1 = self.get_magnitude(color.magnitudes[0], observed=observed)
        m2, s_m2 = self.get_magnitude(color.magnitudes[1], observed=observed)

        # TODO: correlated errors?
        if s_m1 is not None and s_m2 is not None:
            s = np.sqrt(s_m1**2 + s_m2**2)

        return m1 - m2, s
