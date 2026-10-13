PROPOSALID = 'S26B-OT02'
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006
CATID_SCIENCE_CO = 10091
CATID_SCIENCE_GA = 10092
CATID_SCIENCE_GE = 10093

extra_columns = {
    'proposalid': dict(
        constant = PROPOSALID,
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_S26B_{field}_{target_list}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    # Override the minimum number of calibration targets
    netflow_options = dict(
        target_classes = {
            'sky': dict(
                prefix = 'sky',
                min_targets = 400,
                max_targets = 420,
            ),
            'cal': dict(
                prefix = 'cal',
                min_targets = 60,
                max_targets = 240,
            ),
            # Very likely members on the RGB + DEIMOS stars
            "sci_P0": {
                "prefix": "sci",
                # "non_observation_cost": 1000,
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.7
            "sci_P1": {
                "prefix": "sci",
                # "non_observation_cost": 775,
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.0
            "sci_P2": {
                "prefix": "sci",
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0
            },
            "sci_P2.1": {
                "prefix": "sci",
                "non_observation_cost": 900,
                "partial_observation_cost": 100000.0
            },
            "sci_P2.2": {
                "prefix": "sci",
                "non_observation_cost": 800,
                "partial_observation_cost": 100000.0
            },
            "sci_P2.3": {
                "prefix": "sci",
                "non_observation_cost": 700,
                "partial_observation_cost": 100000.0
            },
            "sci_P2.4": {
                "prefix": "sci",
                "non_observation_cost": 600,
                "partial_observation_cost": 100000.0
            },
            # BHB + AGB + ToRGB
            "sci_P3": {
                "prefix": "sci",
                "non_observation_cost": 465,
                "partial_observation_cost": 100000.0
            },
            # Blue Stragglers
            "sci_P4": {
                "prefix": "sci",
                "non_observation_cost": 360,
                "partial_observation_cost": 100000.0
            },
            # Anc 0, Pace, Sestito
            "sci_P5": {
                "prefix": "sci",
                "non_observation_cost": 600,
                "partial_observation_cost": 100000.0
            },
            # DART
            "sci_P5.1": {
                "prefix": "sci",
                "non_observation_cost": 300,
                "partial_observation_cost": 100000.0
            },
            # APOGEE
            "sci_P5.2": {
                "prefix": "sci",
                "non_observation_cost": 200,
                "partial_observation_cost": 100000.0
            },
            # Jensen24
            "sci_P5.3": {
                "prefix": "sci",
                "non_observation_cost": 150,
                "partial_observation_cost": 100000.0
            },
            # Anc 1
            "sci_P6": {
                "prefix": "sci",
                "non_observation_cost": 450,
                "partial_observation_cost": 100000.0
            },
            # Anc 2
            "sci_P7": {
                "prefix": "sci",
                "non_observation_cost": 350,
                "partial_observation_cost": 100000.0
            },
            # Faint dSph members
            "sci_P8": {
                "prefix": "sci",
                "non_observation_cost": 130,
                "partial_observation_cost": 100000.0
            },
            # Fillers
            "sci_P9": {
                "prefix": "sci",
                "non_observation_cost": 10,
                "partial_observation_cost": 100000.0
            }
        },
        cobra_groups = {
            'cal_location': dict(
                min_targets = 3,
            )
        }
    ),

    targets = {
        "hsc": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/sample/SSP/sculptor_SSP_003/hsc_scl_priorities.feather",
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
                    "i_hsc": dict(
                        mag = 'obs_hsc_i',
                        mag_err = 'err_hsc_i',
                    ),
                    "nb515_hsc": dict(
                        mag = 'obs_hsc_nb515',
                        mag_err = 'err_hsc_nb515',
                    ),
                }
            ),
        ),
        "anc": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/scl.anc.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            value_map = {
                'priority': {
                    0: 5,
                    1: 6,
                    2: 7
                }
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_ps1": dict(
                        mag = 'obs_ps_g',
                        mag_err = 'err_ps_g',
                    ),
                    "r_ps1": dict(
                        mag = 'obs_ps_r',
                        mag_err = 'err_ps_r',
                    ),
                    "i_ps1": dict(
                        mag = 'obs_ps_i',
                        mag_err = 'err_ps_i',
                    ),
                    "z_ps1": dict(
                        mag = 'obs_ps_z',
                        mag_err = 'err_ps_z',
                    ),
                },
                limits = {
                    'ps1_g': [16, 23],
                    'ps1_i': [16, 23],
                }
            )
        ),
        "sestito23": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/Sculptor_Sestito23.csv",
            reader_args = dict(),
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 0,
                    dtype = 'int',
                ),
                'exp_time': dict(
                    lambda_args = ['G'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            column_map = {
                'source_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'pmra': 'pmra',
                'pmdec': 'pmdec',
                'EBV': 'ebv',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
                'RV': 'rv',
            },
            photometry = dict(
                filters = {
                    'g_gaia': dict(
                        mag = 'G',
                    ),
                    'bp_gaia': dict(
                        mag = 'BP',
                    ),
                    'rp_gaia': dict(
                        mag = 'RP',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),        
        "dart": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/Sculptor_DART.csv",
            reader_args = dict(),
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5.1,
                    dtype = 'float',
                ),
                'exp_time': dict(
                    lambda_args = ['G'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            column_map = {
                'source_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'EBV': 'ebv',
                'vel': 'rv',
            },
            photometry = dict(
                filters = {
                    'g_gaia': dict(
                        mag = 'G',
                    ),
                    'bp_gaia': dict(
                        mag = 'BP',
                    ),
                    'rp_gaia': dict(
                        mag = 'RP',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),        
        "apogee": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/Sculptor_APOGEE_minusDART.csv",
            reader_args = dict(),
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5.2,
                    dtype = 'float',
                ),
                'exp_time': dict(
                    lambda_args = ['G'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            column_map = {
                'source_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'EBV': 'ebv',
                'pmra': 'pmra',
                'pmdec': 'pmdec',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
                'VHELIO_AVG': 'rv',
            },
            photometry = dict(
                filters = {
                    'g_gaia': dict(
                        mag = 'G',
                    ),
                    'bp_gaia': dict(
                        mag = 'BP',
                    ),
                    'rp_gaia': dict(
                        mag = 'RP',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),        
        "jensen24": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/Sculptor_J24candidates_minus_DART_APOGEE_Sestito23.csv",
            reader_args = dict(),
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5.3,
                    dtype = 'float',
                ),
                'exp_time': dict(
                    lambda_args = ['phot_g_mean_mag'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            column_map = {
                'source_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'EBV': 'ebv',
                'pmra': 'pmra',
                'pmdec': 'pmdec',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
            },
            photometry = dict(
                filters = {
                    'g_gaia': dict(
                        mag = 'phot_g_mean_mag',
                    ),
                    'bp_gaia': dict(
                        mag = 'phot_bp_mean_mag',
                    ),
                    'rp_gaia': dict(
                        mag = 'phot_rp_mean_mag',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),        
        # Miho
        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/sky_sculptor.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = CATID_SKY_PS1,
            extra_columns = extra_columns,
        ),
        # Miho
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sculptor/fluxstd_sculptor.feather",
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax_error': 'err_parallax',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
            },
            prefix = "cal",
            frame = 'icrs',
            epoch = 2016.0,
            catid = CATID_FLUXSTD,
            extra_columns = extra_columns,
            photometry = dict(
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                     # Column storing filter names
                        # psf_mag = f'psf_mag_{b}',                 # All None in data file
                        # psf_mag_err = f'psf_mag_error_{b}',       # All None in data file
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'grizy'
                },
                limits = {
                    'gaia_g': [15, 20.5],
                    'gaia_bp-gaia_g': [0.15, 0.35],
                }
            )
        ),
    }
)
