from typing import List, Dict

from ..config import Config
from .fieldconfig import FieldConfig
from .pointingconfig import PointingConfig
from .targetlistconfig import TargetListConfig
from .netflowoptionsconfig import NetflowOptionsConfig
from ..instrument.instrumentoptionsconfig import InstrumentOptionsConfig
from .gurobioptionsconfig import GurobiOptionsConfig
from .exportoptionsconfig import ExportOptionsConfig
from .debugoptionsconfig import DebugOptionsConfig

class NetflowConfig(Config):
    def __init__(self,
                 field: FieldConfig = FieldConfig(),
                 pointings: List[PointingConfig] = None,
                 targets: Dict[str, TargetListConfig] = {},
                 netflow_options: NetflowOptionsConfig = NetflowOptionsConfig(),
                 instrument_options: InstrumentOptionsConfig = InstrumentOptionsConfig(),
                 gurobi_options: GurobiOptionsConfig = GurobiOptionsConfig(),
                 export_options: ExportOptionsConfig = ExportOptionsConfig(),
                 debug_options: DebugOptionsConfig = DebugOptionsConfig()):
        
        self.field = field
        self.pointings = pointings
        self.targets = targets
        self.netflow_options = netflow_options
        self.instrument_options = instrument_options
        self.gurobi_options = gurobi_options
        self.export_options = export_options
        self.debug_options = debug_options
        
        super().__init__()

    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        config = NetflowConfig(
            netflow_options = NetflowOptionsConfig.default(),
            instrument_options = InstrumentOptionsConfig.default(),
            gurobi_options = GurobiOptionsConfig.default(),
            export_options = ExportOptionsConfig.default(),
            debug_options = DebugOptionsConfig.default(),
        )
        
        return config