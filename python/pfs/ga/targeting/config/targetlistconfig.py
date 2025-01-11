from .config import Config

class TargetListConfig(Config):
    def __init__(self):
        self.path = None
        self.reader = None
        self.reader_args = None
        self.columns = None
        self.column_map = None
        self.data_types = None
        self.index = None
        self.prefix = None
        self.epoch = None
        self.catid = None
        self.proposalid_pattern = None
        self.obcode_pattern = None
        self.filters = None
        self.bands = None
        self.limits = None

        super().__init__()