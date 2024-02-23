from .observation import Observation

class TargetList(Observation):
    def __init__(self, name=None, orig=None):
        super().__init__(name=name, orig=orig)