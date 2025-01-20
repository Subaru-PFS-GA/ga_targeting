import os
from datetime import datetime, timedelta
import numpy as np

import pfs.ga.targeting

config = dict(
    field = dict(
        key = "cosmology_region3_3h",
        name = "cosmology_region3_3h test",
        arms = ['b', 'm', 'n'],
        nvisits = 1,
        exp_time = 30 * 30, # sec
        obs_time = datetime(2025, 8, 20, 11, 0, 0),
        resolution = 'm'
    ),
    pointings = [
        dict(ra=330.682912, dec=-0.518103, posang=30),
    ],
    targets = {
        "targets": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/CO/cosmology_region3_3h_targets.ecsv',
            reader_args = dict(
                comment='#',
                delimiter='\s+'
            ),
            column_map = {
                # 'ID': 'obcode',
                # 'obj_id': 'targetid',
                'R.A.': 'RA',
                'Dec.': 'Dec',
                'Exposure Time': 'exp_time',
                'Priority': 'priority',
            },
            prefix = "sci",
            # epoch = "J2000.0",
            catid = 19002,
            extra_columns = {
                'targetid': dict(
                    lambda_func = R'lambda ID: int(ID[5:])',
                    lambda_args = 'ID',
                    dtype = int,
                ),
                'proposalid': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            },
            photometry = dict(
                filters = {
                    # "g_gaia": dict(
                    #     flux = 'g_gaia'
                    # ),
                    # "g_ps1": dict(
                    #     flux = 'g_ps1'
                    # ),
                    # "g_hsc": dict(
                    #     mag = 'gmag',
                    # ),
                    # "r_hsc": dict(
                    #     mag = 'rmag',
                    # ),
                }
            ),
        ),
        "sky": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/CO/sky_region3_3h_targets.ecsv',
            reader_args = dict(
                comment='#',
                delimiter='\s+'
            ),
            column_map = {
                # 'sky_id': 'targetid',
                # 'ID': 'obcode',
                'R.A.': 'RA',
                'Dec.': 'Dec'
            },
            prefix = "sky",
            extra_columns = {
                'targetid': dict(
                    lambda_func = R'lambda ID: int(float(ID[5:]))',
                    lambda_args = 'ID',
                    dtype = int,
                ),
                'proposalid': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            }
        ),
        "fluxstd": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/CO/star_region3_3h_targets.ecsv',
            reader_args = dict(
                comment='#',
                delimiter='\s+'
            ),
            column_map = {
                # 'ID': 'obcode',
                'R.A.': 'RA',
                'Dec.': 'Dec',
                'Exposure Time': 'exp_time',
            },
            prefix = "cal",
            extra_columns = {
                'targetid': dict(
                    lambda_func = R'lambda ID: int(ID[6:])',
                    lambda_args = 'ID',
                    dtype = int,
                ),
                'proposalid': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_CO_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            },
            photometry = dict(
                filters = {
                    "g_ps1": dict(
                        mag = 'gPS1',
                        flux = "gFluxJy",
                        flux_err = "gFluxJy_err",
                    ),
                    "r_ps1": dict(
                        mag = 'rPS1',
                        flux = "rFluxJy",
                        flux_err = "rFluxJy_err",
                    ),
                    "i_ps1": dict(
                        mag = 'iPS1',
                        flux = "iFluxJy",
                        flux_err = "iFluxJy_err",
                    ),
                    "z_ps1": dict(
                        mag = 'zPS1',
                        flux = "zFluxJy",
                        flux_err = "zFluxJy_err",
                    ),
                    "y_ps1": dict(
                        mag = 'yPS1',
                        flux = "yFluxJy",
                        flux_err = "yFluxJy_err",
                    ),
                }
            ),
        ),
    },
    instrument_options = dict(
        layout = 'calibration',
        cobra_coach_dir = '/tmp/cobra_coach',
        # cobra_coach_module_version = None,
        # instdata_path = None,
        # blackdots_path = None,
        # fiberids_path = None,
        # black_dot_radius_margin = 1.0,
        # spectrograph_modules = [1, 2, 3, 4],
    ),
    netflow_options = dict(
        black_dot_penalty = R'lambda dist: 10 / (dist + 1)',
        cobra_move_cost = R'lambda dist: 5 * dist',
        collision_distance = 2.0,
        forbidden_targets = [],
        forbidden_pairs = [],
        target_classes = {
            'sky': dict(
                prefix = 'sky',
                min_targets = 240,
                max_targets = 320,
                non_observation_cost = 0,
            ),
            'cal': dict(
                prefix = 'cal',
                min_targets = 40,
                max_targets = 240,
                non_observation_cost = 0,
            ),
            'sci_P1': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 1000,
            ),
        },
        cobra_groups = {
            'cal_location': dict(
                # groups = np.random.randint(4, size=2394),
                target_classes = [ 'cal' ],
                min_targets = 3,
                max_targets = 60,
                non_observation_cost = 1000,
            ),
            'sky_instrument': dict(
                # groups = np.random.randint(8, size=2394),
                target_classes = [ 'sky' ],
                min_targets = 10,
                max_targets = 60,
                non_observation_cost = 100,
            ),
            'sky_location': dict(
                # groups = np.random.randint(8, size=2394),
                target_classes = [ 'sky' ],
                min_targets = 10,
                max_targets = 60,
                non_observation_cost = 100,
            )
        },
        time_budgets = None,
        fiber_non_allocation_cost = 1e5,
        num_reserved_fibers = 0,
        allow_more_visits = True,
        epoch = 2016,
        use_named_variables = True,
    ),
    gurobi_options = dict(
        # seed = 0,
        presolve = -1,
        method = 3,
        degenmoves = 0,
        heuristics = 0.5,
        mipfocus = 1,           
        mipgap = 0.01,
        LogToConsole = 1,
        timelimit = 300 # sec
    ),
    debug_options = dict(
        ignore_endpoint_collisions = False,
        ignore_elbow_collisions = False,
        ignore_broken_cobra_collisions = False,
        ignore_forbidden_targets = False,
        ignore_forbidden_pairs = False,
        ignore_calib_target_class_minimum = False,
        ignore_calib_target_class_maximum = False,
        ignore_science_target_class_minimum = False,
        ignore_science_target_class_maximum = False,
        ignore_time_budget = False,
        ignore_cobra_group_minimum = False,
        ignore_cobra_group_maximum = False,
        ignore_reserved_fibers = False,
        ignore_proper_motion = True,
    ),
)
