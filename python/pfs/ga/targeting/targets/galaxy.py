import astropy.units as u

from ..util.args import *
from ..instrument import SubaruHSC
from ..photometry import Photometry, Magnitude, Color
from ..probabilitymap import ProbabilityMap
from .target import Target

class Galaxy(Target):
    def __init__(self,
                 ID,
                 pos,
                 rad=None,
                 pointings=None,
                 **kwargs):
        
        super().__init__(ID, pos, **kwargs)

        # Bounding radius
        self.__rad = normalize_angle(rad, u.arcmin)

        # Telescope pointings with position and position angle
        self.__pointings = pointings

    #region Properties

    def __get_rad(self):
        return self.__rad
    
    rad = property(__get_rad)

    #endregion

    def get_pointings(self, instrument):
        return self.__pointings[instrument]

    def get_photometry(self, instrument):
        raise NotImplementedError()

    def get_cmd(self, instrument):
        raise NotImplementedError()
    
    def get_ccd(self, instrument):
        raise NotImplementedError()
    
    def _get_hsc_dered_mags_colors(self, catalog, mask):
        hsc = SubaruHSC.photometry()
        [ (g0, _), (i0, _), (gi0, _) ] = catalog.get_diagram_values([
                    hsc.magnitudes['g'],
                    hsc.magnitudes['i'],
                    Color([hsc.magnitudes['g'], hsc.magnitudes['i']])
                ], observed=True, mask=mask)
        
        return g0, i0, gi0
    
    def assign_probabilities(self, catalog, pmap, population_id=-1, mask=None):
         # Membership probability
        lp_member, lp_member_mask = pmap.lookup_lp_member(catalog, mask=mask)

        ix = np.where(mask)[0]

        catalog.data['p_member'] = np.nan
        catalog.data['p_member'].iloc[ix[lp_member_mask]] = np.exp(lp_member[:, population_id][lp_member_mask])

    def assign_priorities(self, catalog, mask=None):
        raise NotImplementedError()