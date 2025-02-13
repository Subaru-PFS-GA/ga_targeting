from ..config import Config

class CobraGroupConfig(Config):
    def __init__(self,
                 groups=None,
                 target_classes=None,
                 min_targets=0,
                 max_targets=None,
                 non_observation_cost=10):
        
        # The list of group labels for each cobra, indexed by cobraId - 1
        self.groups = groups

        # The list of target classes that this group is associated with
        self.target_classes = target_classes

        # The minimum number of targets in each cobra group
        self.min_targets = min_targets

        # The maximum number of targets in each cobra group
        self.max_targets = max_targets

        # The cost of not observing a target in this group
        self.non_observation_cost = non_observation_cost

        super().__init__()