import os
from datetime import datetime, timedelta
import numpy as np

import pfs.ga.targeting

path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b25_ENG/ga_targets_outerdisk_l180_b25_ENG_faint-v2.csv'

column_map = {
    'ob_code': 'obcode',
    'obj_id': 'targetid',
    'ra': 'RA',
    'dec': 'Dec',
    'exptime': 'exp_time',
}

extra_columns = {
    'proposalid': dict(
        pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    field = dict(
        key = "outerdisk_l180_b25_bright",
        name = "Outer Disk l=180 b=25 Bright",
        obs_time = datetime(2025, 1, 25, 10, 0, 0),
        exp_time = 60 * 60, # sec
    ),
    pointings = [
        dict(ra=114.4, dec=39.2, posang=30),
    ],
    targets = {
        "gaia": dict(
            path = path,
            mask = 'lambda df: df["input_catalogs"].isin(["SEGUE", "Gaia", "J-PLUS", "Pristine-Gaia", "LAMOST"])',
            column_map = column_map,
            prefix = "sci",
            # epoch = "J2000.0",
            catid = 10088,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        flux = 'g_gaia'
                    )
                }
            )
        ),
        "ps1": dict(
            path = path,
            mask = 'lambda df: df["input_catalogs"] == "PS1"',
            column_map = column_map,
            prefix = "sci",
            # epoch = "J2000.0",
            catid = 10088,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_ps1": dict(
                        flux = 'g_ps1'
                    )
                }
            )
        ),
        "sky": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b25_ENG/outerdisk_b25_sky.feather',
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec'
            },
            prefix = "sky",
            catid = 2007,
            extra_columns = extra_columns,
        ),
        "fluxstd": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b25_ENG/outerdisk_b25_fluxstd.feather',
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'obj_id': 'orig_objid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax': 'parallax',
                'parallax_error': 'err_parallax',
                'pmra': 'pmra',
                'pmra_error': 'err_pmra',
                'pmdec': 'pmdec',
                'pmdec_error': 'err_pmdec',
            },
            mask = 'lambda df: df["prob_f_star"] > 0.5',
            prefix = "cal",
            catid = 3006,
            extra_columns = extra_columns,
            photometry = dict(
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                   # Column storing filter names
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'grizy'
                },
                limits = {
                    'ps1_g': [14, 17],
                }
            )
        ),
    },
)
