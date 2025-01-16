import os
from datetime import datetime, timedelta
import numpy as np

import pfs.ga.targeting

config = dict(
    field = dict(
        key = "mwod_l180_b22_bright",
        name = "Outer Disk l=180 b=22 Bright",
        arms = ['b', 'm', 'n'],
        nvisits = 1,
        exp_time = 30 * 60, # sec
        obs_time = datetime(2025, 1, 25, 10, 0, 0),
        resolution = 'm'
    ),
    pointings = [
        dict(ra=111, dec=38.4, posang=30),
    ],
    targets = {
        "bright": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b22/ga_targets_outerdisk_l180_b22_ENG_bright.csv',
            reader_args = dict(),
            column_map = {
                'ob_code': 'obcode',
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'exptime': 'exp_time',
            },
            prefix = "sci",
            # epoch = "J2000.0",
            catid = 15002,
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
            filters = {
                "g_gaia": dict(
                    flux = 'g_gaia'
                ),
                "g_ps1": dict(
                    flux = 'g_ps1'
                ),
                "g_hsc": dict(
                    mag = 'gmag',
                ),
                "r_hsc": dict(
                    mag = 'rmag',
                ),
            }
        ),
        'vmp': dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/VMP/jplus_field_111_38.csv',
            reader_args = dict(),
            column_map = {
                'SOURCEID': 'targetid',
                'RA': 'RA',
                'DEC': 'Dec',
                'exptime': 'exp_time',
                'PARALLAX': 'parallax',
                'PMRA': 'pmra',
                'PMDEC': 'pmdec',
            },
            prefix = "sci",
            catid = 15003,
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
            filters = {
                "g_gaia": dict(
                    mag = 'GMAG',
                    mag_err = 'GERR'
                ),
                "bp_gaia": dict(
                    mag = 'BPMAG',
                    mag_err = 'BPERR'
                ),
                "rp_gaia": dict(
                    mag = 'RPMAG',
                    mag_err = 'RPERR'
                ),
            }
        ),
        "sky": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b22/outerdisk_b22_sky.feather',
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec'
            },
            prefix = "sky",
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
        ),
        "fluxstd": dict(
            path = '$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_l180_b22/outerdisk_b22_fluxstd.feather',
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
            prefix = "cal",
            proposalid_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
            obcode_pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
            bands = {
                b: dict(
                    filter = f'filter_{b}',                   # Column storing filter names
                    psf_flux = f'psf_flux_{b}',
                    psf_flux_err = f'psf_flux_error_{b}',
                ) for b in 'grizy'
            },

            limits = {
                'ps1_g': [13, 15],
            }
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
        black_dot_penalty = lambda dist: 10 / (dist + 1),
        cobra_move_cost = lambda dist: 5 * dist,
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
            'sci_P0': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 1000,
            ),
            'sci_P3': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 100,
            ),
            'sci_P5': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 10,
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
        LogToConsole = 0,
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
