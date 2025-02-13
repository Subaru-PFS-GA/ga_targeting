from ..config import Config

class TimeBudgetConfig(Config):
    def __init__(self,
                 target_classes=None,
                 budget=5.0):
        
        # List of the target classes this time budget applies to
        self.target_classes = target_classes

        # Total obsevation time of targets in this class, expressed in hours
        self.budget = budget

        super().__init__()