from typing import List, Dict

from pfs.ga.common.config import Config

from .fieldconfig import FieldConfig
from .pointingconfig import PointingConfig
from .targetlistconfig import TargetListConfig
from .netflowoptionsconfig import NetflowOptionsConfig
from ..instrument.instrumentoptionsconfig import InstrumentOptionsConfig
from .gurobioptionsconfig import GurobiOptionsConfig
from .debugoptionsconfig import DebugOptionsConfig
from ..export.exportoptionsconfig import ExportOptionsConfig

class NetflowConfig(Config):
    def __init__(self,
                 field: FieldConfig = FieldConfig(),
                 pointings: List[PointingConfig] = None,
                 targets: Dict[str, TargetListConfig] = {},
                 netflow_options: NetflowOptionsConfig = NetflowOptionsConfig(),
                 instrument_options: InstrumentOptionsConfig = InstrumentOptionsConfig(),
                 gurobi_options: GurobiOptionsConfig = GurobiOptionsConfig(),
                 debug_options: DebugOptionsConfig = DebugOptionsConfig(),
                 export_options: ExportOptionsConfig = ExportOptionsConfig()
    ):
        
        self.field = field
        self.pointings = pointings
        self.targets = targets
        self.netflow_options = netflow_options
        self.instrument_options = instrument_options
        self.gurobi_options = gurobi_options
        self.debug_options = debug_options
        self.export_options = export_options
        
        super().__init__()

    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        config = NetflowConfig()        
        return config