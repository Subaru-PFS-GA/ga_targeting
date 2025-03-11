from typing import Dict
import numpy as np

import pfs.utils
from pfs.utils.fibers import fiberHoleFromFiberId

from ...instrument import SubaruHSC, SubaruPFI
from ..config import Config, Lambda
from .targetclassconfig import TargetClassConfig
from .timebudgetconfig import TimeBudgetConfig
from .cobragroupconfig import CobraGroupConfig

class NetflowOptionsConfig(Config):

    SKY_MIN_TARGETS = 240
    SKY_MAX_TARGETS = 320
    FLUXSTD_MIN_TARGETS = 40
    FLUXSTD_MAX_TARGETS = 240

    def __init__(self,
                 black_dot_penalty: Lambda = None,
                 cobra_move_cost: Lambda = None,
                 target_classes: Dict[str, TargetClassConfig] = None,
                 cobra_groups: Dict[str, CobraGroupConfig] = None,
                 time_budgets: Dict[str, TimeBudgetConfig] = None):
        
        # Add a penalty if the target is too close to a black dot
        self.black_dot_penalty = black_dot_penalty
        # self.black_dot_penalty = lambda dist: 0

        # Do not penalize cobra moves with respect to cobra center
        self.cobra_move_cost = cobra_move_cost
        # self.cobra_move_cost = lambda dist: 0

        self.collision_distance = 2.0
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

        # Convert all source catalogs to this epoch when proper motions are provided
        self.epoch = None       

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

        # TODO: figure out how to make this callable
        #       from object-specific config files

        from ..instrument.instrumentoptionsconfig import InstrumentOptionsConfig

        instrument_options = InstrumentOptionsConfig.default()
        pfi = SubaruPFI(instrument_options=instrument_options)

        cobra_location_labels = pfi.generate_cobra_location_labels(ntheta=6)
        cobra_instrument_labels = pfi.generate_cobra_instrument_labels(ngroups=8)     # for each spectrograph

        cobra_location_labels_count = cobra_location_labels.max() + 1
        cobra_instrument_labels_count = cobra_instrument_labels.max() + 1

        cobra_groups = {
            'sky_location': CobraGroupConfig(
                groups = cobra_location_labels,
                target_classes = [ 'sky' ],
                min_targets = int(np.floor(self.SKY_MIN_TARGETS / cobra_location_labels_count)),
                max_targets = int(np.ceil(self.SKY_MAX_TARGETS / cobra_location_labels_count)),
                non_observation_cost = 100,
            ),
            'sky_instrument': CobraGroupConfig(
                groups = cobra_instrument_labels,
                target_classes = [ 'sky' ],
                min_targets = int(np.floor(self.SKY_MIN_TARGETS / cobra_instrument_labels_count)),
                max_targets = int(np.ceil(self.SKY_MAX_TARGETS / cobra_instrument_labels_count)),
                non_observation_cost = 100,
            ),
            'cal_location': CobraGroupConfig(
                groups = cobra_location_labels,
                target_classes = [ 'cal' ],
                min_targets = int(np.floor(self.FLUXSTD_MIN_TARGETS / cobra_location_labels_count)),
                max_targets = int(np.ceil(self.FLUXSTD_MAX_TARGETS / cobra_location_labels_count)),
                non_observation_cost = 1000,
            ),
            # 'cal_instrument': CobraGroupConfig(
            #     groups = cobra_instrument_labels,
            #     target_classes = [ 'cal' ],
            #     min_targets = int(np.floor(self.FLUXSTD_MIN_TARGETS / cobra_instrument_labels_count)),
            #     max_targets = int(np.ceil(self.FLUXSTD_MAX_TARGETS / cobra_instrument_labels_count)),
            #     non_observation_cost = 1000,
            # ),
        }

        return cobra_groups

    def __create_target_classes(self, num_classes=10):
        target_classes = {
            'sky': TargetClassConfig(
                prefix = 'sky',
                min_targets = self.SKY_MIN_TARGETS,
                max_targets = self.SKY_MAX_TARGETS,
                non_observation_cost = 0,
            ),
            'cal': TargetClassConfig(
                prefix = 'cal',
                min_targets = self.FLUXSTD_MIN_TARGETS,
                max_targets = self.FLUXSTD_MAX_TARGETS,
                non_observation_cost = 0,
            ),
        }
            
        non_observation_cost = np.ceil(np.logspace(2, 3, num_classes)).astype(int)[::-1]
        for i in range(num_classes):
            target_classes[f'sci_P{i}'] = TargetClassConfig(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = non_observation_cost[i],
                partial_observation_cost = 1e5,
            )

        return target_classes