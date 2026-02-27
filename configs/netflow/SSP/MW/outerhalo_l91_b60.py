from datetime import datetime


DATA_DIR = '$PFS_TARGETING_DATA/data/targeting/MW/outerhalo_l91_b60_SSP'

PROPOSALID = 'S26A-OT02'
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3011
CATID_SCIENCE_CO = 10091
CATID_SCIENCE_GA = 10092
CATID_SCIENCE_GE = 10093

ID_PREFIX = 0x5E000000000

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
        pattern = "SSP_GA_{obs_time:%Y%m%d}_{field}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    field = dict(
        key = "outerhalo_l91_b60_faint",
        name = "GA Outer Halo l=91 b=60 Faint",
        obs_time = datetime(2026, 3, 17, 6, 0, 0),
        id_prefix = ID_PREFIX
    ),
    pointings = [
        dict(ra=217.6569, dec=50.4532, posang=0.0, priority=1),
    ],
    netflow_options = dict(
       cobra_groups = {
           'cal_location': dict(
               min_targets = 0,
           ),
       }
    ),
    targets = {
        # Miho
        "ps1": dict(
            path = f'{DATA_DIR}/ga_targets_outerhalo_l91_b60_faint.csv',
            mask = 'lambda df: df["input_catalogs"] == "PS1"',
            column_map = column_map,
            value_map = {
                'priority': {
                    1: 1,
                    2: 2,
                    3: 3,
                    4: 4,
                    5: 5,
                    8: 8,
                    9: 9
                }
            },
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
                #limits = {
                #    'ps1_g-ps1_r': [-0.5, 1.5],
                #}
            )
        ),
        "sdssv": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/CC/SDSS-V/minesweeper_v1.0.0_mag16.0.feather",
            #mask = 'lambda df: df["input_catalogs"] == "Gaia"',
            column_map = column_map,
            value_map = {
                'priority': {
                    1: 0
                }
            },
            prefix = "sci",
            # epoch = "J2000.0",
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "bp_gaia": dict(
                        mag = 'BP',
                    ),
                    "g_gaia": dict(
                        mag = 'G',
                    ),
                    "rp_gaia": dict(
                        mag = 'RP',
                    )
                },
                bands = {
                    b: dict(
                        filter = f'filter_{b}',
                        psf_flux = f'psf_flux_{b}',
                        #psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'gri'
                },
                #limits = {
                #    'ps1_g-ps1_r': [-0.5, 1.5],
                #}
            )
        ),
        "segue": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/CC/SEGUE/segue_sspparam_mag16.5.feather",
            #mask = 'lambda df: df["input_catalogs"] == "Gaia"',
            column_map = column_map,
            value_map = {
                'priority': {
                    1: 0
                }
            },
            prefix = "sci",
            # epoch = "J2000.0",
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "bp_gaia": dict(
                        mag = 'phot_bp_mean_mag',
                    ),
                    "g_gaia": dict(
                        mag = 'phot_g_mean_mag',
                    ),
                    "rp_gaia": dict(
                        mag = 'phot_rp_mean_mag',
                    )
                },
                bands = {
                    b: dict(
                        filter = f'filter_{b}',
                        psf_flux = f'psf_flux_{b}',
                        #psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'gri'
                },
                #limits = {
                #    'ps1_g-ps1_r': [-0.5, 1.5],
                #}
            )
        ),
        "sky": dict(
            path = f'{DATA_DIR}/l91b60_sky.csv',
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
            path = f'{DATA_DIR}/l91b60_fluxstd.csv',
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
            extra_columns = {
                **{
                    'exp_time': dict(
                        constant = 0.0,
                        dtype = 'float'
                    ),
                    'priority': dict(
                        constant = -1,
                        dtype = 'int'
                    )
                },
                **extra_columns
            },
            photometry = dict(
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                   # Column storing filter names
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'grizy'
                },
                limits = {
                    'ps1_g': [16.5, 20.0],
                }
            )
        ),
    },
)

