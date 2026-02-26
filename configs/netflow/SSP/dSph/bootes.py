PROPOSALID = 'S25A-OT02'
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
                min_targets = 0,
            )
        }
    ),
    targets = {
        "hsc": dict(
            path = "$PFS_TARGETING_DATA/dSph/bootes/sample/SSP/bootes/hsc_bootes_priorities.feather",
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
                    "g_sdss": dict(
                        mag = 'obs_sdss_g',
                        mag_err = 'err_sdss_g',
                    ),
                    "r_sdss": dict(
                        mag = 'obs_sdss_r',
                        mag_err = 'err_sdss_r',
                    ),
                #     "g_hsc": dict(
                #         mag = 'obs_hsc_g',
                #         mag_err = 'err_hsc_g',
                #     ),
                #     "i_hsc": dict(
                #         mag = 'obs_hsc_i',
                #         mag_err = 'err_hsc_i',
                #     ),
                #     "nb515_hsc": dict(
                #         mag = 'obs_hsc_nb515',
                #         mag_err = 'err_hsc_nb515',
                #     ),
                }
            )
        ),
        "sky": dict(
            path = "$PFS_TARGETING_DATA/dSph/bootes/sky_bootesi_v2.feather",
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

        # MIHO NEW
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/dSph/bootes/fluxstd_bootesi_v2.feather",
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax_error': 'err_parallax',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
            },
            mask = 'lambda df: df["prob_f_star"] > 0.5',
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
                    'ps1_g': [16, 19],
                    'ps1_g-ps1_r': [0.25, 0.4],
                }
            )
        ),
    },
)
