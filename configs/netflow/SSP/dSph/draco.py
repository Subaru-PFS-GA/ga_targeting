import numpy as np
from datetime import datetime, timedelta

from pfs.ga.targeting.instrument import SubaruHSC

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
        pattern = "PFS_SSP_GA_S25A_{name}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    debug_options = dict(
        ignore_calib_target_class_minimum = False,
    ),
    netflow_options = dict(
        target_classes = {
            'sky': dict(
                prefix = 'sky',
                min_targets = 240,
                max_targets = 320,
            ),
            'cal': dict(
                prefix = 'cal',
                min_targets = 40,
                max_targets = 240,
            ),
            # Very likely members on the RGB + DEIMOS stars
            "sci_P0": {
                "prefix": "sci",
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.7
            "sci_P1": {
                "prefix": "sci",
                "non_observation_cost": 775,
                "partial_observation_cost": 100000.0
            },
            # Bright p_member > 0.0
            "sci_P2": {
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
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/sample/draco_nb/hsc_dra_priorities.feather",
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
                    "i2_hsc": dict(
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
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/dra.anc-short_exposures.feather",
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
        "jpass_vmp": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/draco_JPLUS_VMP.csv",
            column_map =
            {
                'SOURCEID': 'targetid',
                'RA': 'RA',
                'DEC': 'Dec',
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5,
                    dtype = 'int',
                ),
                'exp_time': dict(
                    lambda_args = ['GMAG'],
                    lambda_func = "lambda g0: 1800 * np.maximum(np.minimum(np.rint(5 * ((g0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'GMAG',
                    ),
                    "bp_gaia": dict(
                        mag = 'BPMAG',
                    ),
                    "rp_gaia": dict(
                        mag = 'RPMAG',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),
        "jpass_more": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/draco_JPLUS_more.csv",
            column_map =
            {
                'SOURCEID': 'targetid',
                'RA': 'RA',
                'DEC': 'Dec',
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 6,
                    dtype = 'int',
                ),
                'exp_time': dict(
                    lambda_args = ['GMAG'],
                    lambda_func = "lambda g0: 1800 * np.maximum(np.minimum(np.rint(5 * ((g0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'GMAG',
                    ),
                    "bp_gaia": dict(
                        mag = 'BPMAG',
                    ),
                    "rp_gaia": dict(
                        mag = 'RPMAG',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),
        "jingkun": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/dra_jingkun.csv",
            column_map =
            {
                'name': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5,
                    dtype = 'int',
                ),
                'exp_time': dict(
                    lambda_args = ['phot_g_mean_mag'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'phot_g_mean_mag',
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),
        "pristine": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/Dra_VMP_G16.fits",
            column_map = {
                'source_id': 'targetid',
                'Plx': 'parallax',
                'pmRA': 'pmra',
                'pmDE': 'pmdec',
            },
            prefix = "sci",
            frame= 'icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = {
                **extra_columns,
                'priority': dict(
                    constant = 5,
                    dtype = 'int',
                ),
                'exp_time': dict(
                    lambda_args = ['RPmag'],
                    lambda_func = "lambda r0: 1800 * np.maximum(np.minimum(np.rint(5 * ((r0 - 16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)",
                    dtype = 'int'
                )
            },
            photometry = dict(
                filters = {
                    "g_gaia": dict(
                        mag = 'Gmag',
                        mag_err = 'e_Gmag'
                    ),
                    "bp_gaia": dict(
                        mag = 'BPmag',
                        mag_err = 'e_BPmag'
                    ),
                    "rp_gaia": dict(
                        mag = 'RPmag',
                        mag_err = 'e_RPmag'
                    ),
                },
                limits = {
                    'gaia_g': [16, 23],
                }
            )
        ),
        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/sky_draco.feather",
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

        # # DOBOS PS1
        # "fluxstd": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/PS1_Dra_fluxstd_3_dobos.feather",
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
        #     catid = CATID_SCIENCE_GA,
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
        #         },
        #         limits = {
        #             'ps1_g': [16, 19],
        #             'ps1_i': [16, 19],
        #             'ps1_g-ps1_r': [0.25, 0.4],
        #         }
        #     )
        # ),

        # MIHO NEW
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/fluxstd_draco_miho_20250315.feather",
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
        #     path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/guide_draco.feather",
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
