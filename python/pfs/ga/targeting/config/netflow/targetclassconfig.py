from pfs.ga.common.config import Config

class TargetClassConfig(Config):
    def __init__(self,
                 prefix=None,
                 min_targets=0,
                 max_targets=None,
                 non_observation_cost=10,
                 partial_observation_cost=1e5):
        
        self.prefix = prefix
        self.min_targets = min_targets
        self.max_targets = max_targets
        self.non_observation_cost = non_observation_cost
        self.partial_observation_cost = partial_observation_cost

        super().__init__()