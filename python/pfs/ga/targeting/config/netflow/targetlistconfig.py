from typing import Dict

from pfs.ga.common.config import Config, Lambda

from .extracolumnconfig import ExtraColumnConfig
from .photometryconfig import PhotometryConfig

class TargetListConfig(Config):
    def __init__(self,
                 extra_columns: Dict[str, ExtraColumnConfig] = None,
                 photometry: PhotometryConfig = None,
                 mask: Lambda = None):
        
        # Path to the data file
        self.path = None

        # Data reader class to use
        self.reader = None

        # Arguments to pass to the data reader
        self.reader_args = None

        # List of columns in the file, use when columns cannot be inferred or to override
        self.columns = None

        # Mapping of columns to their names, key is original name, value is new name
        self.column_map = None

        # Data types for each column, if cannot be inferred. Otherwise, override
        # the data types. Passed to `DataFrame.astype`
        self.data_types = None

        # Dict of dict to map values in columns. Passed to `DataFrame.replace`.
        self.value_map = None

        # A dict of extra column definition to be added to the dataset after loading.
        self.extra_columns = extra_columns

        # Index column or list of index columns. Passed to `DataFrame.set_index`
        self.index = None

        # Dataset prefix, 'sci', 'cal' or 'sky'
        self.prefix = None

        # Catalog coordinate frame, defaults to ICRS
        self.frame = None

        # Catalog coordinate system equinox, default to None when using ICRS
        # otherwise must be specified
        self.equinox = None

        # Overrides the epoch column of the datasets. Sets epoch of coordinates, use for
        # high proper motion stars. This is not the equinox of the coordinates.
        self.epoch = None

        # Overrides that catId that will be written to the design file.
        self.catid = None

        # Overrides the priority of the targets
        self.priority = None

        # Overrides the exposure time of the targets
        self.exp_time = None

        # Definition of the filters for the photometric system
        self.photometry = photometry

        # Catalog filter, expressed as a lambda, rendered in string
        self.mask = None

        # Radius in arcseconds to cross-match sources with the rest of the catalogs
        self.crossmatch_radius = None

        super().__init__()