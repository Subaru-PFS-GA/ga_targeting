import os
from typing import List, Dict
import numpy as np

import pfs.utils
from pfs.utils.fibers import fiberHoleFromFiberId

from .config import Config
from .fieldconfig import FieldConfig
from .targetlistconfig import TargetListConfig
from ..instrument import SubaruHSC, SubaruPFI

class NetflowConfig(Config):
    def __init__(self,
                 field: FieldConfig = FieldConfig(),
                 targets: Dict[str, TargetListConfig] = {}):
        
        self.field = field
        self.targets = targets

        self.instrument_options = dict(
            # Instrument layout mode
            layout = 'calibration',     # Use calibration product
            # layout = 'full',          # Use default settings

            # Temp directory for cobra coach output
            cobra_coach_dir = '/tmp/cobra_coach',

            # Cobra coach module version
            # cobra_coach_module_version = None,

            # FPI configuration
            # instdata_path = None,

            # Black dot list
            # blackdots_path = None,

            # Black dot radius margin
            # black_dot_radius_margin = 1.0,

            # List of spectrograph modules to use
            # spectrograph_modules = [1, 2, 3, 4],
        )
        
        self.netflow_options = dict(
            # Add a penalty if the target is too close to a black dot
            black_dot_penalty = None,
            # black_dot_penalty = lambda dist: 0,

            fiber_non_allocation_cost = 1e5,

            collision_distance = 2.0,
            elbow_collisions = True,
            # forbidden_targets = [
            #     43218108431
            # ],
            # forbidden_pairs = [
            #     [43486543901, 43218108431],
            # ],

            target_classes = self.__create_target_classes(),
            cobra_groups = self.__create_cobra_groups(self.instrument_options),

            # time_budgets = {
            #     'science': dict(
            #         target_classes = [ 'sci_P0', 'sci_P1', 'sci_p2' ],
            #         budget = 5  # hr
            #     )
            # },

            # Do not penalize cobra moves with respect to cobra center
            cobra_move_cost = lambda dist: 0,

            num_reserved_fibers = 0,

            # This will only be used when netflow is rewritten to step-by-step fiber assignment
            # constrain_already_observed = False,

            # Allow more visits than minimally required
            allow_more_visits = True,

            epoch = 2016, # all catalogs must match
            ignore_proper_motion = True,

            # FPI configuration
            fiberids_path = os.path.join(os.path.dirname(pfs.utils.__file__), '../../../data/fiberids'),

            # Generate full gurobi variable names instead of numbered ones (slow to build problem)
            use_named_variables = False,
        )

        self.gurobi_options = dict(
            seed=0,                 # random seed
            presolve=-1,            # agressiveness of presolve which tries to eliminate variables from the LP problem
            method=3,               # 3 means concurrent, 4 means deterministic concurrent
            degenmoves=0,           # degenerate simplex moves, set to 0 to prevent too much time to be spent on trying to improve the current solution
            heuristics=0.5,         # how much of the time to spend by performing heuristics
            mipfocus=1,             # mipfocus=1 is balanced toward finding more feasible solutions
                                    # mipfocus=2 is balanced toward proving that the current solution is the best
                                    # mipfocus=3 is to be used when the objection bound is moving very slowly
            mipgap=0.01,            # relative stopping criterion for bounds on the objective
            LogToConsole=1,         # 
            timelimit=300           # in sec
        )

        self.debug_options = dict(
            ignore_endpoint_collisions = False,
            ignore_elbow_collisions = False,
            ignore_forbidden_pairs = False,
            ignore_forbidden_singles = False,
            ignore_calib_target_class_minimum = False,
            ignore_calib_target_class_maximum = False,
            ignore_science_target_class_minimum = False,
            ignore_science_target_class_maximum = False,
            ignore_time_budget = False,
            ignore_cobra_group_minimum = False,
            ignore_cobra_group_maximum = False,
            ignore_reserved_fibers = False,
        )

        super().__init__()

    def __create_cobra_groups(self, instrument_options):
        """
        Generate the cobra group definitions to distribute the sky fiber uniformly
        across the focal plane as well as the slit.
        """

        pfi = SubaruPFI(instrument_options=instrument_options)

        cobra_location_labels = self.__create_cobra_location_labels(pfi, ntheta=6)
        cobra_instrument_labels = self.__create_cobra_instrument_labels(pfi, ngroups=16)

        cobra_groups = {
            'sky_location': dict(
                # groups = np.random.randint(8, size=ncobras),
                groups = cobra_location_labels,
                target_classes = [ 'sky' ],
                min_targets = 20,
                max_targets = 80,
                non_observation_cost = 100,
            ),
            'sky_instrument': dict(
                # groups = np.random.randint(8, size=ncobras),
                groups = cobra_instrument_labels,
                target_classes = [ 'sky' ],
                min_targets = 3,
                max_targets = 20,
                non_observation_cost = 100,
            ),
            'cal_location': dict(
                # groups = np.random.randint(8, size=ncobras),
                groups = cobra_location_labels,
                target_classes = [ 'cal' ],
                min_targets = 5,
                max_targets = 10,
                non_observation_cost = 1000,
            ),
            # 'cal_instrument': dict(
            #     # groups = np.random.randint(8, size=ncobras),
            #     groups = cobra_instrument_labels,
            #     target_classes = [ 'cal' ],
            #     min_targets = 10,
            #     max_targets = 60,
            #     non_observation_cost = 1000,
            # ),
        }

        return cobra_groups
    
    def __create_cobra_location_labels(self, pfi, ntheta=6):

        # Get the focal plane coordinates of the cobras
        x, y = pfi.bench.cobras.centers[:].real, pfi.bench.cobras.centers[:].imag

        # Convert to polar coordinates around the center of the focal plane
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)

        # Assign labels to the cobras based on the polar coordinate
        theta_bins = np.linspace(-np.pi, np.pi, ntheta + 1)
        r_bins = np.array([0, 150, 240])

        theta_labels = np.digitize(theta, theta_bins, right=False) - 1
        r_labels = np.digitize(r, r_bins, right=False) - 1
        cobra_location_labels = (r_bins.size - 1) * theta_labels + r_labels

        # Add one more label in the center
        cobra_location_labels[r < 60] = cobra_location_labels.max() + 1

        return cobra_location_labels
    
    def __create_cobra_instrument_labels(self, pfi, ngroups=8):
        ncobras = len(pfi.bench.cobras.centers)
        cobra_instrument_labels = np.zeros(ncobras, dtype=int)
        mask = pfi.fiber_map.cobraId != 65535

        fiber_hole_group_size = pfi.fiber_map.fiberHoleId.max() / ngroups
        cobra_instrument_labels[pfi.fiber_map.cobraId[mask] - 1] = \
            (pfi.fiber_map.spectrographId[mask] - 1) * ngroups \
            + (np.round(pfi.fiber_map.fiberHoleId[mask] - 1) / fiber_hole_group_size).astype(int)

        return cobra_instrument_labels

    def __create_target_classes(self):
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
                calib = True,
            ),
        }
            
        for i in range(10):
            target_classes[f'sci_P{i}'] = dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = max(20 - 2 * i, 1),
                partial_observation_cost = 1e5,
            )

        target_classes[f'sci_P0']['non_observation_cost'] = 1000
        target_classes[f'sci_P1']['non_observation_cost'] = 500
        target_classes[f'sci_P2']['non_observation_cost'] = 200
        target_classes[f'sci_P3']['non_observation_cost'] = 100
        target_classes[f'sci_P4']['non_observation_cost'] = 100
        target_classes[f'sci_P5']['non_observation_cost'] = 100
        target_classes[f'sci_P6']['non_observation_cost'] = 100
        target_classes[f'sci_P7']['non_observation_cost'] = 50
        target_classes[f'sci_P8']['non_observation_cost'] = 10
        target_classes[f'sci_P9']['non_observation_cost'] = 0

        return target_classes