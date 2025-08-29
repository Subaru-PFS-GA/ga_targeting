from pfs.ga.common.config import Config

class InstrumentOptionsConfig(Config):
    def __init__(self):
        
        # Use the default values ('full') or load the calibration data for the PFI ('calibration')
        self.layout = 'calibration'

        # Temp directory for cobra coach output
        self.cobra_coach_dir = '/tmp/cobra_coach'

        # Version of the cobra coach module to be used
        self.cobra_coach_module_version = None

        # Path to the instrument data, defaults to the one under the `pfs.instdata` package
        self.instdata_path = None

        # Black dot list, defaults to the one under the `pfs.instdata` package
        self.blackdots_path = None

        # FPI configuration, defaults to the one under the `pfs.instdata` package
        self.fiberids_path = None

        # The margin factor in radius of the black dots to avoid in fiber allocation.
        # This parameter is passed to the `Bench` object.
        self.black_dot_radius_margin = 1.0

        # List of the spectgraph modules to be used.
        self.spectrograph_modules = [1, 2, 3, 4]

        # Ignore any errors in the calibration products
        self.ignore_calibration_errors = False

        super().__init__()