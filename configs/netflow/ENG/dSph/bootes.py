from datetime import datetime

PROPOSALID = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}"
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006
CATID_SCIENCE_GA = 10088

extra_columns = {
    'proposalid': dict(
        constant = PROPOSALID,
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{field}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    field = dict(
        key = "bootes",
        name = "GA Bootes I dSph ENG",
        obs_time = datetime(2025, 1, 25, 10, 0, 0),
        exp_time = 6 * 30 * 60, # sec
    ),
    pointings = [
        dict(ra=210.025, dec=14.5, posang=30, exp_time=6 * 30 * 60),
    ],
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
                min_targets = 0,
                max_targets = 20,
            )
        },
        epoch = 2016.0,
    ),
    targets = {
        "dsph": dict(
            path = f"$PFS_TARGETING_DATA/data/targeting/dSph/bootesi/bootesi_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            frame='icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
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
            path = f"$PFS_TARGETING_DATA/data/targeting/dSph/bootesi/sky_bootesi.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = CATID_SKY_PS1,
            extra_columns = extra_columns,
        ),
        "fluxstd": dict(
            path = f"$PFS_TARGETING_DATA/data/targeting/dSph/bootesi/fluxstd_bootesi.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "cal",
            frame = 'icrs',
            epoch = 2016.0,
            catid = CATID_FLUXSTD,
            extra_columns = extra_columns,
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
)
