FIELD = 'PFS_22'
PROPOSALID = 'S25B-OT02'
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
        pattern = "SSP_GA_S25B_{sector}_{target_list}_{{targetid:d}}_{resolution}",
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
            # PNe
            "sci_P0": {
                "prefix": "sci",
                # "non_observation_cost": 1000,
                "non_observation_cost": 2500,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.95
            "sci_P1": {
                "prefix": "sci",
                # "non_observation_cost": 775,
                "non_observation_cost": 1500,
                "partial_observation_cost": 100000.0
            },
            # Faint p_member > 0.95
            "sci_P2": {
                "prefix": "sci",
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.80
            "sci_P3": {
                "prefix": "sci",
                "non_observation_cost": 465,
                "partial_observation_cost": 100000.0
            },
            # Faint p_member > 0.80
            "sci_P4": {
                "prefix": "sci",
                "non_observation_cost": 360,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.1
            "sci_P5": {
                "prefix": "sci",
                "non_observation_cost": 600,
                "partial_observation_cost": 100000.0
            },
            # Faint p_member > 0.1
            "sci_P6": {
                "prefix": "sci",
                "non_observation_cost": 450,
                "partial_observation_cost": 100000.0
            },
            # Anc 1
            "sci_P7": {
                "prefix": "sci",
                "non_observation_cost": 350,
                "partial_observation_cost": 100000.0
            },
            # Anc 2
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
            path = "$PFS_TARGETING_DATA/data/targeting/m31/m31_all_SSP/sample/m31_all_SSP/hsc_m31_priorities.feather",
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
            )
        ),
        "anc": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/m31/m31.ucd.anc.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            value_map = {
                'priority': {
                    0: 7,
                    1: 7,
                    2: 8
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
        "pne": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/m31/m31.pne.feather",
            # reader = None
            reader_args = dict(),
            #column_map = {'objid': 'targetid'},
            value_map = {
                'priority': {
                    0: 0,
                }
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_cfht": dict(
                        mag = 'g_cfht',
                    ),
                },
                limits = {
                    'cfht_g': [16, 24.5],
                }
            )
        ),

        # # Use the PS1 x GAIA sample as a fall-back for empty fibers
        # # This entry also add PS1 and GAIA photometry to all other objects
        # "gaia": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/m31/PS1_GAIA_point_source_M31_2_dobos.csv",
        #     column_map = {
        #         'source_id': 'targetid',
        #         'ra': 'RA',
        #         'dec': 'Dec',
        #         'radial_velocity': 'rv',
        #         'ref_epoch': 'epoch',
        #     },
        #     prefix = "sci",
        #     frame= 'icrs',
        #     epoch = 2016.0,
        #     catid = CATID_SCIENCE_GA,
        #     extra_columns = {
        #         **extra_columns,
        #         'priority': dict(
        #             constant = 9,
        #             dtype = 'int',
        #         ),
        #         'exp_time': dict(
        #             lambda_args = ['rPSFMag'],
        #             lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
        #             dtype = 'int'
        #         )
        #     },
        #     photometry = dict(
        #         filters = {
        #             "g_gaia": dict(
        #                 flux = 'phot_g_mean_flux',
        #                 flux_err = 'phot_g_mean_flux_error'
        #             ),
        #             "bp_gaia": dict(
        #                 flux = 'phot_bp_mean_flux',
        #                 flux_err = 'phot_bp_mean_flux_error'
        #             ),
        #             "rp_gaia": dict(
        #                 flux = 'phot_rp_mean_flux',
        #                 flux_err = 'phot_rp_mean_flux_error'
        #             ),
        #             "g_ps1": dict(
        #                 mag = 'gPSFMag',
        #                 mag_err = 'gPSFMagErr',
        #             ),
        #             "r_ps1": dict(
        #                 mag = 'rPSFMag',
        #                 mag_err = 'rPSFMagErr',
        #             ),
        #             "i_ps1": dict(
        #                 mag = 'iPSFMag',
        #                 mag_err = 'iPSFMagErr',
        #             ),
        #             "z_ps1": dict(
        #                 mag = 'zPSFMag',
        #                 mag_err = 'zPSFMagErr',
        #             ),
        #             "y_ps1": dict(
        #                 mag = 'yPSFMag',
        #                 mag_err = 'yPSFMagErr',
        #             ),
        #         },
        #         limits = {
        #             'gaia_rp': [16, 23],
        #         }
        #     )
        # ),

        "sky": dict(
            path = f"$PFS_TARGETING_DATA/data/targeting/m31/M31_fluxstd_sky/m31_sky.feather",
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

        # DOBOS
        # "fluxstd": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/PS1_UMI_fluxstd_2_dobos.feather",
        #     reader_args = dict(),
        #     column_map = {
        #         'obj_id': 'targetid',
        #         # 'ra': 'RA',
        #         # 'dec': 'Dec',
        #         'parallax_error': 'err_parallax',
        #         'pmra_error': 'err_pmra',
        #         'pmdec_error': 'err_pmdec',
        #         'radial_velocity': 'rv',
        #         'radial_velocity_error': 'err_rv',
        #     },
        #     # mask = 'lambda df: df["prob_f_star"] > 0.5',
        #     # mask = 'lambda df: df["prob_f_star"] > 0.1',
        #     prefix = "cal",
        #     frame = 'icrs',
        #     epoch = 2016.0,
        #     catid = SCIENCE_GA_CATID,
        #     extra_columns = extra_columns,
        #     photometry = dict(
        #         filters = {
        #             'g_ps1': dict(
        #                 mag = 'gPSFMag',
        #                 mag_err = 'gPSFMagErr'
        #             ),
        #             'r_ps1': dict(
        #                 mag = 'rPSFMag',
        #                 mag_err = 'rPSFMagErr'
        #             ),
        #             'i_ps1': dict(
        #                 mag = 'iPSFMag',
        #                 mag_err = 'iPSFMagErr'
        #             ),
        #             'z_ps1': dict(
        #                 mag = 'zPSFMag',
        #                 mag_err = 'zPSFMagErr'
        #             ),
        #             'y_ps1': dict(
        #                 mag = 'yPSFMag',
        #                 mag_err = 'yPSFMagErr'
        #             )
        #         }
        #     )
        # ),

        # MIHO NEW
        "fluxstd": dict(
            path = f"$PFS_TARGETING_DATA/data/targeting/m31/M31_fluxstd_sky/m31_fluxstd.feather",
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

        # "guide": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/m31/guide_m31.feather",
        #     reader_args = dict(),
        #     prefix = "ag",
        #     photometry = dict(
        #         filters = {
        #             "gaia_g": dict(
        #                 flux = "flux_gaia_g",
        #                 flux_err = "err_flux_gaia_g",
        #             ),
        #             "gaia_bp": dict(
        #                 flux = "flux_gaia_bp",
        #                 flux_err = "err_flux_gaia_bp",
        #             ),
        #             "gaia_rp": dict(
        #                 flux = "flux_gaia_rp",
        #                 flux_err = "err_flux_gaia_rp",
        #             )
        #         }
        #     )
        # )
    },
)
