from typing import List, Dict

from pfs.ga.common.config import Config

from .fieldconfig import FieldConfig
from .pointingconfig import PointingConfig
from .targetlistconfig import TargetListConfig
from .netflowoptionsconfig import NetflowOptionsConfig
from ..instrument.instrumentoptionsconfig import InstrumentOptionsConfig
from .gurobioptionsconfig import GurobiOptionsConfig
from .debugoptionsconfig import DebugOptionsConfig

class NetflowConfig(Config):
    def __init__(self,
                 field: FieldConfig = FieldConfig(),
                 pointings: List[PointingConfig] = None,
                 targets: Dict[str, TargetListConfig] = {},
                 netflow_options: NetflowOptionsConfig = NetflowOptionsConfig(),
                 instrument_options: InstrumentOptionsConfig = InstrumentOptionsConfig(),
                 gurobi_options: GurobiOptionsConfig = GurobiOptionsConfig(),
                 debug_options: DebugOptionsConfig = DebugOptionsConfig()):
        
        self.field = field
        self.pointings = pointings
        self.targets = targets
        self.netflow_options = netflow_options
        self.instrument_options = instrument_options
        self.gurobi_options = gurobi_options
        self.debug_options = debug_options
        
        super().__init__()

    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        config = NetflowConfig()
        config.netflow_options = NetflowOptionsConfig.default()
        config.instrument_options = InstrumentOptionsConfig.default()
        config.gurobi_options = GurobiOptionsConfig.default()
        config.debug_options = DebugOptionsConfig.default()
        
        return config