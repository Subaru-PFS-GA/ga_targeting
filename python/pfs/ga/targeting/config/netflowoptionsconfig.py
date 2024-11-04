from typing import Dict
import numpy as np

import pfs.utils
from pfs.utils.fibers import fiberHoleFromFiberId

from ..instrument import SubaruHSC, SubaruPFI
from .config import Config
from .targetclassconfig import TargetClassConfig
from .timebudgetconfig import TimeBudgetConfig
from .cobragroupconfig import CobraGroupConfig

class NetflowOptionsConfig(Config):
    def __init__(self,
                 target_classes: Dict[str, TargetClassConfig] = None,
                 cobra_groups: Dict[str, CobraGroupConfig] = None,
                 time_budgets: Dict[str, TimeBudgetConfig] = None):
        
        # Add a penalty if the target is too close to a black dot
        self.black_dot_penalty = None
        # self.black_dot_penalty = lambda dist: 0

        # Do not penalize cobra moves with respect to cobra center
        self.cobra_move_cost = None
        # self.cobra_move_cost = lambda dist: 0

        self.collision_distance = 2.0
        self.elbow_collisions = True
        self.forbidden_targets = []
        self.forbidden_pairs = [
            # [43486543901, 43218108431],
        ]

        self.target_classes = target_classes
        self.cobra_groups = cobra_groups
        self.time_budgets = time_budgets

        # Penalty for not allocating a fiber to a target
        self.fiber_non_allocation_cost = 1e5

        # Reserve fibers
        self.num_reserved_fibers = 0

        # Allow more visits than minimally required
        self.allow_more_visits = True

        self.epoch = 2016 # all catalogs must match

        # Generate full gurobi variable names instead of numbered ones (slow to build problem)
        # It is necessary only, when the netflow solution will be loaded from a file to
        # extract assignments based on variable names.
        self.use_named_variables = True

        super().__init__()
    
    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        config = NetflowOptionsConfig()
        config.target_classes = config.__create_target_classes()
        config.cobra_groups = config.__create_cobra_groups()
        return config

    def __create_cobra_groups(self):
        """
        Generate the cobra group definitions to distribute the sky fiber uniformly
        across the focal plane as well as the slit.
        """

        from .instrumentoptionsconfig import InstrumentOptionsConfig

        instrument_options = InstrumentOptionsConfig.default()
        pfi = SubaruPFI(instrument_options=instrument_options)

        cobra_location_labels = self.__create_cobra_location_labels(pfi, ntheta=6)
        cobra_instrument_labels = self.__create_cobra_instrument_labels(pfi, ngroups=16)

        cobra_groups = {
            'sky_location': CobraGroupConfig(
                groups = cobra_location_labels,
                target_classes = [ 'sky' ],
                min_targets = 20,
                max_targets = 80,
                non_observation_cost = 100,
            ),
            'sky_instrument': CobraGroupConfig(
                groups = cobra_instrument_labels,
                target_classes = [ 'sky' ],
                min_targets = 3,
                max_targets = 20,
                non_observation_cost = 100,
            ),
            'cal_location': CobraGroupConfig(
                groups = cobra_location_labels,
                target_classes = [ 'cal' ],
                min_targets = 3,
                max_targets = 20,
                non_observation_cost = 1000,
            ),
            # 'cal_instrument': CobraGroupConfig(
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
            'sky': TargetClassConfig(
                prefix = 'sky',
                min_targets = 240,
                max_targets = 320,
                non_observation_cost = 0,
            ),
            'cal': TargetClassConfig(
                prefix = 'cal',
                min_targets = 40,
                max_targets = 240,
                non_observation_cost = 0,
            ),
        }
            
        for i in range(10):
            target_classes[f'sci_P{i}'] = TargetClassConfig(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = max(20 - 2 * i, 1),
                partial_observation_cost = 1e5,
            )

        target_classes[f'sci_P0'].non_observation_cost = 1000
        target_classes[f'sci_P1'].non_observation_cost = 500
        target_classes[f'sci_P2'].non_observation_cost = 200
        target_classes[f'sci_P3'].non_observation_cost = 100
        target_classes[f'sci_P4'].non_observation_cost = 100
        target_classes[f'sci_P5'].non_observation_cost = 100
        target_classes[f'sci_P6'].non_observation_cost = 100
        target_classes[f'sci_P7'].non_observation_cost = 50
        target_classes[f'sci_P8'].non_observation_cost = 10
        target_classes[f'sci_P9'].non_observation_cost = 0

        return target_classes