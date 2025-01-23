from datetime import datetime, timedelta
from pfs.ga.targeting.targets.dsph import Fornax
from pfs.ga.targeting.instrument import SubaruHSC

DATA_DIR = '/datascope/subaru'

config = dict(
    field = dict(
        key = "bootes",
        name = "GA Bootes I dSph ENG",
        obs_time = datetime(2025, 1, 25, 10, 0, 0),
        exp_time = 6 * 30 * 60, # sec
    ),
    pointings = [
        dict(ra=210.025, dec=14.5, posang=30),
    ],
    targets = {
        "dsph": dict(
            path = f"{DATA_DIR}/data/targeting/dSph/bootesi/bootesi_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            epoch = "J2000.0",
            catid = 10088,
            extra_columns = {
                'proposalid': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            },
            photometry = dict(
                filters = {
                    "g_hsc": dict(
                        mag = 'obs_hsc_g',
                        mag_err = 'err_hsc_g',
                    ),
                    "r_hsc": dict(
                        mag = 'obs_hsc_r',
                        mag_err = 'err_hsc_r',
                    ),
                }
            ),
        ),
        "sky": dict(
            path = f"{DATA_DIR}/data/targeting/dSph/bootesi/sky_bootesi.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = 1007,
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
            catid = 3006,
            extra_columns = {
                'proposalid': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            },
            photometry = dict(
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
