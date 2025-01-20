from datetime import datetime, timedelta
from pfs.ga.targeting.targets.dsph import Fornax
from pfs.ga.targeting.instrument import SubaruHSC

DATA_DIR = '/datascope/subaru'

config = dict(
    targets = {
        "dsph": dict(
            path = f"{DATA_DIR}/data/targeting/dSph/bootesi/bootesi_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            epoch = "J2000.0",
            catid = 15001,
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
            filters = {
                "g_hsc": dict(
                    mag = 'obs_hsc_g',
                    mag_err = 'err_hsc_g',
                    # flux = "g_hsc_flux",
                    # flux_err = "g_hsc_flux_err",
                    # psf_mag = "g_hsc",
                    # psf_mag_err = "g_hsc_err",
                    # psf_flux = "g_hsc_flux",
                    # psf_flux_err = "g_hsc_flux_err",
                    # fiber_mag = "g_hsc",
                    # fiber_mag_err = "g_hsc_err",
                    # fiber_flux = "g_hsc_flux",
                    # fiber_flux_err = "g_hsc_flux_err",
                    # total_mag = "g_hsc",
                    # total_mag_err = "g_hsc_err",
                    # total_flux = "g_hsc_flux",
                    # total_flux_err = "g_hsc_flux_err",
                ),
                "r_hsc": dict(
                    mag = 'obs_hsc_r',
                    mag_err = 'err_hsc_r',
                ),
            }
        ),
        "sky": dict(
            path = f"{DATA_DIR}/data/targeting/dSph/bootesi/sky_bootesi.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
        ),
        "fluxstd": dict(
            path = f"{DATA_DIR}/data/targeting/dSph/bootesi/fluxstd_bootesi.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "cal",
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
            bands = {
                b: dict(
                    filter = f'filter_{b}',                   # Column storing filter names
                    psf_mag = f'psf_mag_{b}',
                    psf_mag_err = f'psf_mag_error_{b}',
                    psf_flux = f'psf_flux_{b}',
                    psf_flux_err = f'psf_flux_error_{b}',
                ) for b in 'gri'
            },
            limits = {
                'ps1_g': [17, 19],
            }
        ),
    },
    # Override the minimum number of calibration targets
    netflow_options = dict(
        target_classes = {
            'cal': dict(
                min_targets = 200,
                max_targets = 300,
            )
        },
        cobra_groups = {
            'cal_location': dict(
                min_targets = 15,
                max_targets = 20,
            )
        }
    )
)
