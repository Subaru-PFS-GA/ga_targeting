from .config import Config

class DebugOptionsConfig(Config):
    def __init__(self):
        self.ignore_endpoint_collisions = False
        self.ignore_elbow_collisions = False
        self.ignore_broken_cobra_collisions = False
        self.ignore_forbidden_targets = False
        self.ignore_forbidden_pairs = False
        self.ignore_calib_target_class_minimum = False
        self.ignore_calib_target_class_maximum = False
        self.ignore_science_target_class_minimum = False
        self.ignore_science_target_class_maximum = False
        self.ignore_time_budget = False
        self.ignore_cobra_group_minimum = False
        self.ignore_cobra_group_maximum = False
        self.ignore_reserved_fibers = False
        self.ignore_proper_motion = False
        self.ignore_missing_priority = True
        self.ignore_missing_exp_time = True

        super().__init__()