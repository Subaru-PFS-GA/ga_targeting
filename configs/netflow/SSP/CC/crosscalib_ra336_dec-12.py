from datetime import datetime

DATA_DIR = '$PFS_TARGETING_DATA/data/targeting/CC/crosscalib_ra336_dec-12'

# TODO: update these IDs
PROPOSALID = 'S25A-OT02'
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006
CATID_SCIENCE_CO = 10091
CATID_SCIENCE_GA = 10092
CATID_SCIENCE_GE = 10093

ID_PREFIX = 0x40000000000

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
        key = "crosscalib_ra336_decm12",
        name = "Cross-Calibration ra=336 dec=-12",
        obs_time = datetime(2025, 6, 1, 15, 0, 0),
        id_prefix = ID_PREFIX
    ),
    pointings = [
        dict(ra=335.75, dec=-12.00, posang=0.0, priority=200),
    ],
    targets = {
        # Miho
        "gaia": dict(
            path = f'{DATA_DIR}/ga_crosscalib_ra336_dec-12.ecsv',
            # mask = 'lambda df: df["input_catalogs"] == "Gaia"',
            column_map = column_map,
            prefix = "sci",
            # epoch = "J2000.0",
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **{
                    'exp_time': dict(
                        constant = 300,
                        dtype = 'float'
                    ),
                },
                **extra_columns
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        flux = 'g_gaia',
                    ),
                },
                bands = {
                    b: dict(
                        filter = f'filter_{b}',
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'gi'
                },
            )
        ),
        # Jingkun
        "gaiaxp": dict(
            path = f'{DATA_DIR}/gaiaxp_mp_bright.fits',
            column_map = {
                'SOURCE_ID': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sci",
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **{
                    'exp_time': dict(
                        constant = 300,
                        dtype = 'float'
                    ),
                    'priority': dict(
                        constant = 0,
                        dtype = 'int'
                    )
                },
                **extra_columns
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'phot_g_mean_mag',
                    ),
                },
                # limits = {
                #     'gaia_g': [16, 23],
                # }
            )
        ),
        # Jingkun
        "skymapper": dict(
            path = f'{DATA_DIR}/skymapper_mp_bright.fits',
            column_map = {
                'SOURCEID': 'targetid',
                'RA': 'RA',
                'DEC': 'Dec',
            },
            prefix = "sci",
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **{
                    'exp_time': dict(
                        constant = 300,
                        dtype = 'float'
                    ),
                    'priority': dict(
                        constant = 0,
                        dtype = 'int'
                    )
                },
                **extra_columns
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'G_C',
                    ),
                },
                # limits = {
                #     'gaia_g': [16, 23],
                # }
            )
        ),
        # Federico
        "pristine": dict(
            path = f'{DATA_DIR}/Pristine_synth_335.7_-12.0.csv',
            prefix = "sci",
            catid = CATID_SCIENCE_GA,
            priority = 0,
            column_map = {
                'GaiaDR3': 'targetid',
                'RAICRS': 'RA',
                'DEICRS': 'Dec',
            },
            extra_columns = {
                **{
                    'exp_time': dict(
                        constant = 300,
                        dtype = 'float'
                    )
                },
                **extra_columns
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'Gmag0',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ), 
        # Miho
        "sky": dict(
            path = f'{DATA_DIR}/ga_crosscalib_ra336_dec-12_sky.csv',
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
        # Miho
        "fluxstd": dict(
            path = f'{DATA_DIR}/ga_crosscalib_ra336_dec-12_fluxstd.csv',
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
                # limits = {
                #     'ps1_g': [16.5, 18.5],
                # }
            )
        ),
    },
)
