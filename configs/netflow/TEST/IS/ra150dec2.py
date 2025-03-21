import numpy as np
from datetime import datetime, timedelta

from pfs.ga.targeting.instrument import SubaruHSC

PROPOSALID = 'xxx'
CATID_SCIENCE = 9999
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006

extra_columns = {
    'proposalid': dict(
        constant = PROPOSALID,
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "IS_S25A_ra150dec2_{target_list}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    field = dict(
        key = 'ra150dec2',
        id_prefix = 0,
    ),
    debug_options = dict(
        ignore_calib_target_class_minimum = False,
    ),
    pointings = [
        dict(ra=150.05, dec=2.2, posang=0, obs_time='2025-03-25 10:00:00', exp_time=46800, priority=0, stage=0)
    ],
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
            "sci_P0": {
                "prefix": "sci",
                "non_observation_cost": 1000,
                "partial_observation_cost": 100000.0,
            },
            "sci_P1": {
                "prefix": "sci",
                "non_observation_cost": 775,
                "partial_observation_cost": 100000.0,
            },
            "sci_P2": {
                "prefix": "sci",
                "non_observation_cost": 600,
                "partial_observation_cost": 100000.0,
                "max_targets": 450,
            },
            "sci_P3": {
                "prefix": "sci",
                "non_observation_cost": 465,
                "partial_observation_cost": 100000.0,
                "max_targets": 800,
            },
            "sci_P4": {
                "prefix": "sci",
                "non_observation_cost": 360,
                "partial_observation_cost": 100000.0
            },
            "sci_P5": {
                "prefix": "sci",
                "non_observation_cost": 600,
                "partial_observation_cost": 100000.0
            },
            "sci_P6": {
                "prefix": "sci",
                "non_observation_cost": 450,
                "partial_observation_cost": 100000.0
            },
            "sci_P7": {
                "prefix": "sci",
                "non_observation_cost": 350,
                "partial_observation_cost": 100000.0
            },
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
    ),
    targets = {
        "galaxies": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/misc/IS/targets_for_laci.csv",
            # reader = None
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'exptime': 'exp_time'
            },
            prefix = "sci",
            frame='icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    'i2_hsc': dict(
                        flux = 'i2_hsc',
                        flux_err = 'i2_hsc_err'
                    )
                }
            )
        ),

        # # Use the PS1 x GAIA sample as a fall-back for empty fibers
        # # This entry also add PS1 and GAIA photometry to all other objects
        # "gaia": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/dSph/draco/PS1_GAIA_point_source_Draco_2_dobos.csv",
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

        "filler": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/misc/IS/PS1_point_source_ra150dec2_dobos.csv",
            # reader = None
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
            },
            prefix = "sci",
            frame='icrs',
            epoch = 2016.0,
            priority = 0,
            exp_time = 46800,
            catid = CATID_SCIENCE,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    f'{b}_ps1': dict(
                        mag = f'{b}PSFMag',
                        mag_err = f'{b}PSFMagErr'
                    ) for b in 'grizy'
                },
                limits = {
                    'ps1_r' : [20, None]
                }
            )
        ),


        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/misc/IS/ra150dec2_sky.feather",
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

        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/misc/IS/ra150dec2_fluxstd.feather",
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
