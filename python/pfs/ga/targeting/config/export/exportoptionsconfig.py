from typing import List, Dict
from datetime import datetime

from pfs.ga.common.config import Config, Lambda

class ExportOptionsConfig(Config):
    def __init__(
        self,
        export_repeats = None,
        filter_map: Dict[str, str] = None
    ):
        
        # Export each repeat of a field as a separate target list
        # This will create separate ppcCodes but the same pfsDesignId
        self.export_repeats = export_repeats

        # Map from the filter names used in the input catalogs to the filter names
        # used in the output catalogs.
        self.filter_map = filter_map

    @classmethod
    def default(cls):
        return cls(
            export_repeats = True,
            filter_map = {}
        )
