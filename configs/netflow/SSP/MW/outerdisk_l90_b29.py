from datetime import datetime

DATA_DIR = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l90_b29_SSP'

# TODO: update these IDs
PROPOSALID = 'S25A-OT02'
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006
CATID_SCIENCE_CO = 10091
CATID_SCIENCE_GA = 10092
CATID_SCIENCE_GE = 10093

column_map = {
    'ob_code': 'obcode',
    'obj_id': 'targetid',
    'ra': 'RA',
    'dec': 'Dec',
    'exptime': 'exp_time',
}

extra_columns = {
    'proposalid': dict(
        pattern = PROPOSALID,
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{field}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    field = dict(
        key = "outerdisk_l90_b29_faint",
        name = "GA Outer Disk l=90 b=29 Faint",
        obs_time = datetime(2025, 5, 28, 12, 0, 0),
    ),
    pointings = [
        dict(ra=270.725, dec=61.00, posang=120.0, priority=0),
    ],
    netflow_options = dict(
        cobra_groups = {
            'cal_location': dict(
                min_targets = 3,
            ),
        }
    ),
    targets = {
        "ps1": dict(
            path = f'{DATA_DIR}/ga_targets_outerdisk_l90_b29_faint.ecsv',
            mask = 'lambda df: df["input_catalogs"] == "PS1"',
            column_map = column_map,
            prefix = "sci",
            # epoch = "J2000.0",
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_ps1": dict(
                        mag = 'gmag',
                    ),
                    "r_ps1": dict(
                        flux = 'r_ps1',
                        mag = 'rmag'
                    )
                },
                bands = {
                    b: dict(
                        filter = f'filter_{b}',
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'grizy'
                },
            )
        ),
        "sky": dict(
            path = f'{DATA_DIR}/l90b29_sky.csv',
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'catalog_id': 'catid',
                'ra': 'RA',
                'dec': 'Dec'
            },
            prefix = "sky",
            catid = CATID_SKY_PS1,
            extra_columns = extra_columns,
        ),
        "fluxstd": dict(
            path = f'{DATA_DIR}/l90b29_fluxstd.csv',
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'obj_id': 'orig_objid',
                'catalog_id': 'catid',
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
            catid = CATID_FLUXSTD,
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
                    'ps1_g': [16.5, 18.5],
                }
            )
        ),
    },
)
