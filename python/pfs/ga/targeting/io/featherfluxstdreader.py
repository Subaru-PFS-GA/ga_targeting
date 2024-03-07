import pandas as pd

from ..data import Observation
from .observationreader import ObservationReader

class FeatherFluxStdReader(ObservationReader):
    """
    Reads sky coordinates from feather files
    """
    
    def __init__(self, orig=None):
        super().__init__(orig=orig)

        if not isinstance(orig, FeatherFluxStdReader):
            pass
        else:
            pass

    def read(self, filename):
        df = pd.read_feather(filename,)

        df.rename(inplace=True,
                  columns={
                      'fluxstd_id': 'objid',
                      'obj_id': 'orig_objid',
                      'ra': 'RA',
                      'dec': 'Dec',
                      'parallax': 'parallax',
                      'parallax_error': 'err_parallax',
                      'pmra': 'pmra',
                      'pmra_error': 'err_pmra',
                      'pmdec': 'pmdec',
                      'pmdec_error': 'err_pmdec',
                    #   'tract': 'tract',
                    #   'patch': 'patch',
                    #   'target_type_id': 'target_type_id',
                    #   'input_catalog_id': 'input_catalog_id',
                      'psf_mag_g': 'obs_hsc_g',
                      'psf_mag_error_g': 'err_hsc_g',
                      'psf_mag_r': 'obs_hsc_r',
                      'psf_mag_error_r': 'err_hsc_r',
                      'psf_mag_i': 'obs_hsc_i',
                      'psf_mag_error_i': 'err_hsc_i',
                      'psf_mag_z': 'obs_hsc_z',
                      'psf_mag_error_z': 'err_hsc_z',
                      'psf_mag_y': 'obs_hsc_y',
                      'psf_mag_error_y': 'err_hsc_y',
                      'psf_mag_j': 'obs_hsc_j',             # This probably isn't HSC
                      'psf_mag_error_j': 'err_hsc_j',
                      'psf_flux_g': 'obs_flux_hsc_g',
                      'psf_flux_error_g': 'err_flux_hsc_g',
                      'psf_flux_r': 'obs_flux_hsc_r',
                      'psf_flux_error_r': 'err_flux_hsc_r',
                      'psf_flux_i': 'obs_flux_hsc_i',
                      'psf_flux_error_i': 'err_flux_hsc_i',
                      'psf_flux_z': 'obs_flux_hsc_z',
                      'psf_flux_error_z': 'err_flux_hsc_z',
                      'psf_flux_y': 'obs_flux_hsc_y',
                      'psf_flux_error_y': 'err_flux_hsc_y',
                      'psf_flux_j': 'obs_flux_hsc_j',
                      'psf_flux_error_j': 'err_flux_hsc_j',
                    #   'prob_f_star': 'prob_f_star',
                    #   'flags_dist': 'flags_dist',
                    #   'flags_ebv': 'flags_ebv',
                    #   'version': 'version',
                    #   'created_at',
                    #   'updated_at',
                    #   'filter_g',
                    #   'filter_r',
                    #   'filter_i',
                    #   'filter_z',
                    #   'filter_y',
                    #   'filter_j', 
                    #   'teff_brutus',
                    #   'teff_brutus_low',
                    #   'teff_brutus_high',
                    #   'logg_brutus',
                    #   'logg_brutus_low',
                    #   'logg_brutus_high'
                  })
        
        c = self._create_catalog()
        c._set_data(df)

        return c
