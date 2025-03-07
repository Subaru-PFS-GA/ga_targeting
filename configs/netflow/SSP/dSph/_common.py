import os
from datetime import datetime, timedelta
import numpy as np

import pfs.ga.targeting

config = dict(
    field = dict(
        arms = ['b', 'm', 'n'],
        nvisits = 1,
        resolution = 'm'
    ),
    instrument_options = dict(
        layout = 'calibration',
        cobra_coach_dir = '/tmp/cobra_coach',
        # cobra_coach_module_version = None,
        # instdata_path = None,
        # blackdots_path = None,
        # fiberids_path = None,
        black_dot_radius_margin = 1.65,
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
        },
        cobra_groups = {
            'cal_location': dict(
                # groups = np.random.randint(4, size=2394),
                target_classes = [ 'cal' ],
                min_targets = 8,
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
        epoch = 2016.0,
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
        ignore_proper_motion = False,
        ignore_missing_priority = True,
        ignore_missing_exp_time = True,
    ),
)
