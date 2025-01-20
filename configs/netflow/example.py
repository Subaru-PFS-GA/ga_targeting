import os
from datetime import datetime, timedelta
import numpy as np

import pfs.ga.targeting

# This is an example config file that demonstrates all configurable features of
# the netflow algorithm, as implemented by this library. Typically, most settings can
# be left default to successfully run the target selection.
#
# When using the python format, the configuration file works as an executable python
# script, that creates a dictionary in the `config` variable. Alternatively, the
# configuration can be a json or yaml file. When the active configuration is saved
# by the script, it will be saved in the yaml format for readability.
#
# The structure of the configuration dictionary must follow the schema defined in
# the `pfs.ga.targeting.config` namespace starting from the `NetflowConfig` class
# as the root.
#
# We suggest using the `dict()` notation whenever the name of the configuration
# attributes are defined in the schema, as it will be easier to maintain the code.
# When the name of the attributes are not defined in the schema, however, the curly
# braces notation is suggested with string keys.
#
# Note, that when the `ga-netflow` script is executed from the command-line, in addition
# to the configuration file specified in the `--config` argument, a well-known target
# can also be specified with the `--dsph` argument and the script loads many of the
# configuration settings from the python library class of the well-known target. When
# a configuration file is also specified, settings in the configuration take precedence,
# so it is possible to override the boresight of the pointings, for example, with a
# custom configuration file. The `ga-netflow` script can use multiple configuration files,
# which are loaded in the same order as specified on the command-line. Settings in files
# specified later take precedence over settings in files specified earlier.

config = dict(

    # This section defines general properties of the field.
    field = dict(
        # Short name for the main object of the field
        key = "umi",
        
        # Full name of the field
        name = "Ursa Minor Dwarf",

        # List of arms to be used to observe the field. It has no effect on the
        # targeting algorithm, but these values are used to generate the pfsDesign files.
        arms = ['b', 'm', 'n'],

        # Number of visits for each pointing.
        nvisits = 6,

        # Exposure time for each visit, in seconds.
        exp_time = 30 * 60, # sec

        # Obsevation time, UTC, for all visits. Although objects move in the focal plane over
        # time, the difference in position is assumed to be negligible for the purposes of targeting.
        # The pfsDesign files will be further tweaked to position the fibers to the correct positions.
        obs_time = datetime(2025, 5, 1, 10, 0, 0),

        # Specify m for medium resolution or l for low resolution of the red arm
        resolution = 'm'
    ),

    # The following section defines the pointings of the field.
    pointings = [
        dict(ra=227.1, dec=67.25, posang=30),
        dict(ra=227.3, dec=67.25, posang=30),
    ],

    # The following section defines the input target lists, including the science targets,
    # sky fibers, and calibration targets. The configuration must define the path to the
    # target list file, the reader arguments, the column mapping, the target prefix, etc.
    targets = {
        "dsph": dict(
            # Path to the target list file
            path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi.feather'),

            # Reader object to be used to read the target list file. If None, the file will
            # be read by pandas or else, depending on the extension.
            # reader = None

            # Options to be passed to the reader or the pandas read function.
            reader_args = dict(),

            # Column mapping between the columns in the target list file and the columns
            # used by the netflow algorithm. This dictionary is passed to `DataFrame.rename`, so
            # the keys are the current column names and the values are the new column names.
            # The netflow algorithm requires the presence of the following columns:
            # `targetid`, `ra`, `dec`, `epoch`, `priority`, `exp_time`. Optionally, define
            # the following columns: `pmra`, `pmdec`, `parallax`, `class`, `prefix`, `req_visits`,
            # `done_visits`, `penalty`.
            # Column name are case-sensitive.
            column_map = {
                'objid': 'targetid'
            },

            # Replace values in columns
            # value_map = {
            #     'priority': {
            #         1: 2
            #     }
            # }

            # The prefix of the target type. Can be `sci`, `sky`, or `cal`.
            prefix = "sci",

            # Override the epoch of the coordinates in the target list.
            epoch = "J2000.0",

            # Override the catalog ID of the target list
            catid = 15001,

            # Extra columns can be generated by three different methods: i) constant value,
            # ii) f-string pattern, or iii) lambda function. The data type must be specified
            # in all cases with pandas string representation, see:
            # https://pandas.pydata.org/pandas-docs/stable/user_guide/basics.html#basics-dtypes
            # When using lambda expression, the source code must be represented as a string,
            # otherwise it's not seriazable to string when saving the configuration to json or yaml.
            extra_columns = {
                'proposalid': dict(
                    # Override the proposal ID of the target list
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    # Override the observation code pattern of the target list
                    # Single curly braces are resolved first with values that are common
                    # to all objects of the target list. Double curly braces are resolved
                    # second on an object-by-object basis.
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string',
                ),
                # 'constant_column': dict(
                #     # Generate a column with a constant value.
                #     value = 1,
                #     dtype = 'int',
                # ),
                # 'lambda_column': dict(
                #     # Generate a column with a lambda function.
                #     lambda_func = R'lambda x: x + 1',
                #     lambda_args = 'targetid',
                #     dtype = 'int',
                # ),
            },

            # Magnitudes/fluxes and filter names can be provided in two different ways:
            #
            # 1. If the configuration entry `filters` is provided, the filter names are listed in
            #    the configuration file as a dictionary, where the keys are the filter names.
            #    The columns corresponding to the various fluxes and magnitudes must be listed
            #    individually in the configuration.
            #
            # 2. If the configuration entry `bands` is provided, the band names are listed in
            #    the configuration file as a dictionary, where the keys are the band names.
            #    The columns corresponding to the various fluxes and magnitudes must be listed
            #    individually in the configuration for each band. In addition, a column must be
            #    defined that contains the filter name for each band. This method is typically
            #    used for calibration target lists where the targets might come from different
            #    photometric catalogs.
            photometry = dict(
                # Define the filters used to measure the flux of the targets.
                # At least one filter must be defined and either the flux or the magnitude must be
                # defined for each filter. The missing fluxes and/or magnitudes are automatically
                # calculated from the provided values but no aperture or fiber effects are taken
                # into account, so the values are only indicative.
                filters = {
                    "g_hsc": dict(
                        mag = 'sdss_g',
                        mag_err = 'err_sdss_g',
                        # flux = "g_hsc_flux",
                        # flux_err = "g_hsc_flux_err",
                        # psf_mag = "g_hsc",
                        # psf_mag_err = "g_hsc_err",
                        # psf_flux = "g_hsc_flux",
                        # psf_flux_err = "g_hsc_flux_err",
                        # fiber_mag = "g_hsc",
                        # fiber_mag_err = "g_hsc_err",
                        # fiber_flux = "g_hsc_flux",
                        # fiber_flux_err = "g_hsc_flux_err",
                        # total_mag = "g_hsc",
                        # total_mag_err = "g_hsc_err",
                        # total_flux = "g_hsc_flux",
                        # total_flux_err = "g_hsc_flux_err",
                    ),
                    "i_hsc": dict(
                        mag = 'sdss_r',
                        mag_err = 'err_sdss_r',
                    ),
                }
            )
        ),
        "sky": dict(
            path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi_sky.feather'),
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec'
            },
            prefix = "sky",
            extra_columns = {
                'proposalid': dict(
                    # Override the proposal ID of the target list
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    # Override the observation code pattern of the target list
                    # Single curly braces are resolved first with values that are common
                    # to all objects of the target list. Double curly braces are resolved
                    # second on an object-by-object basis.
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string',
                ),
            },
        ),
        "fluxstd": dict(
            path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi_fluxstd.feather'),
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'obj_id': 'orig_objid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax': 'parallax',
                'parallax_error': 'err_parallax',
                'pmra': 'pmra',
                'pmra_error': 'err_pmra',
                'pmdec': 'pmdec',
                'pmdec_error': 'err_pmdec',
            },
            prefix = "cal",
            extra_columns = {
                'proposalid': dict(
                    # Override the proposal ID of the target list
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    # Override the observation code pattern of the target list
                    # Single curly braces are resolved first with values that are common
                    # to all objects of the target list. Double curly braces are resolved
                    # second on an object-by-object basis.
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string',
                ),
            },
            photometry = dict(
                # Define the bands used to measure the flux of the calibration targets. The filter
                # name is taken from the column defined in the `filter` key. Missing fluxes and
                # magnitudes are automatically calculated from the provided values as described above.
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                   # Column storing filter names
                        # psf_mag = f'psf_mag_{b}',
                        # psf_mag_err = f'psf_mag_error_{b}',
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                        # fiber_mag = None,
                        # fiber_mag_err = None,
                        # fiber_flux = None,
                        # fiber_flux_err = None,
                        # total_mag = None,
                        # total_mag_err = None,
                        # total_flux = None,
                        # total_flux_err = None,
                    ) for b in 'gr'
                },

                # Define magnitude limits for the calibration targets. The limits can only include
                # minimum and maximum limits on the magnitudes. If magnitude are not available,
                # the are calculated from the corresponding fluxes. The dictionary keys must match
                # be composed of '{photometry}_{magnitude}'. The names stored in the data frame columns
                # are parsed automatically into the name of the photometric systems and the filter.
                # Magnitude (flux) types are searched in the following order: `psf`, `fiber`, `total`.
                limits = {
                    'ps1_g': [17, 19],
                }
            ),
        ),
    },

    # The following options control the instrument configuration.
    instrument_options = dict(
        # Load the calibration data for the PFI
        layout = 'calibration',
        # Use the default (full) layout for the PFI - for testing only!
        # layout = 'full',

        # Temp directory for cobra coach output
        cobra_coach_dir = '/tmp/cobra_coach',

        # Version of the cobra coach module to be used
        # cobra_coach_module_version = None,

        # Path to the instrument data, defaults to the one under the `pfs.instdata` package
        # instdata_path = None,

        # Black dot list, defaults to the one under the `pfs.instdata` package
        # blackdots_path = None,

        # FPI configuration, defaults to the one under the `pfs.instdata` package
        # fiberids_path = None,

        # The margin factor in radius of the black dots to avoid in fiber allocation.
        # This parameter is passed to the `Bench` object.
        # black_dot_radius_margin = 1.0,

        # List of the spectgraph modules to be used.
        # spectrograph_modules = [1, 2, 3, 4],
    ),

    # The following options control the netflow algorithm.
    netflow_options = dict(
        # Add a penalty if the target is too close to a black dot. It must be None or a
        # lambda expression represented as a string.
        #
        # The typical distance of the targets from a black dot depends on the target density but
        # it is about 3-5 mm for dSph targets such as the Ursa Minor dwarf. The penalty should
        # be set to a value that is comparable to the cost of not observing a science target.
        # Smaller distances should mean larger penalties but penalties must be positive.
        black_dot_penalty = R'lambda dist: 10 / (dist + 1)',

        # Cost of moving a cobra away from its center position. It must be None or a lambda
        # expression represented as a string.
        #
        # The diameter of the patrol region of the cobras is 9.5 mm. Compare this to the
        # typical total cost of the problem and the MIP gap in the solver options.
        # The cost of moving a cobra away from its center position should be comparable,
        # or larger than the remaining gap after optimization, otherwise the cobra move cost
        # will be ignored.
        cobra_move_cost = R'lambda dist: 5 * dist',

        # The minimum distance between cobra tips and cobra elbows to avoid collisions.
        collision_distance = 2.0,

        # List of target ids that are excluded from the target lists.
        forbidden_targets = [],

        # List of pairs (tuples) of target ids that are not allowed to be observed together.
        forbidden_pairs = [],

        # Definitions of target classes. By default, classes are defined for `sky`, `cal` and `sci`.
        # Finer target classes are defined for the well-known targets in the python library.
        # The example given below overrides the default classes. The name of the target class is
        # defined in the dictionary key, and the keys must correspond to the `priority` column in
        # the target lists.
        target_classes = {
            'sky': dict(
                # Primary type of the target. Must be `sky` for sky fibers.
                prefix = 'sky',

                # Minimum number of fibers to be allocated to sky
                min_targets = 240,

                # Maximum number of fibers to be allocated to sky
                max_targets = 320,

                # Cost of not observing a particular sky position. This value should be
                # smaller than the cost of not observing a science target.
                non_observation_cost = 0,
            ),
            'cal': dict(
                # Primary type of the target. Must be `cal` for flux standards.
                prefix = 'cal',
                min_targets = 40,
                max_targets = 240,
                non_observation_cost = 0,
            ),
            'sci_P0': dict(
                # Primary type of the target. Must be `sci` for science target classes.
                prefix = 'sci',

                min_targets = None,
                max_targets = None,

                # The non-observation cost is the primary quantity that controls science target
                # allocation. Set it to a high value for high-priority targets. Define values
                # that increase super-linearly with the priority.
                non_observation_cost = 1000,
            ),
            'sci_P1': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 500,
            ),
            'sci_P2': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 200,
            ),
            'sci_P3': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 100,
            ),
            'sci_P4': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 50,
            ),
            'sci_P5': dict(
                prefix = 'sci',
                min_targets = None,
                max_targets = None,
                non_observation_cost = 10,
            ),
        },

        # Cobra groups define constraint on the minimum and maximum number of fibers in a group
        # that are to be allocated to a specific class of calibration targets, i.e. sky or cal.
        # These settings are used to ensure that the flux standards are more or less uniformly
        # distributed over the field of view and sky fibers are uniformly distributed along the
        # spectrograph slits. Cobra groups are defined by assigning an integer label to each cobra.
        # The example below defines the cobra groups randomly, which is not a good practice. Use
        # the default configuration instead.
        cobra_groups = {
            'cal_location': dict(
                # Cobra group labels, for each cobra
                groups = np.random.randint(4, size=2394),

                # The target prefix this cobra group applies to
                target_classes = [ 'cal' ],

                # The minimum number of fibers in a group that are to be allocated to the target class
                min_targets = 10,

                # The maximum number of fibers in a group that are to be allocated to the target class
                max_targets = 60,

                # The cost of not observing a target in this class
                non_observation_cost = 1000,
            ),
            'sky_instrument': dict(
                groups = np.random.randint(8, size=2394),
                target_classes = [ 'sky' ],
                min_targets = 30,
                max_targets = 40,
                non_observation_cost = 100,
            )
        },

        # Time budgets allocated for a set of target classes.
        time_budgets = {
            'science': dict(
                # List of the target classes this time budget applies to
                target_classes = [ 'sci_P0', 'sci_P1', 'sci_P2', 'sci_P3', 'sci_P4', 'sci_P5' ],

                # Total obsevation time of targets in this class, expressed in hours
                budget = 5 # hr
            )
        },

        # The netflow cost term associated with a fiber that is not allocated to any target.
        fiber_non_allocation_cost = 1e5,

        # The number of reserved fibers that are not allocated to any target.
        num_reserved_fibers = 0,

        # Allow more visits than minimally required. Use this setting to allocate more time
        # to bright targets in less dense fields with multiple visits.
        allow_more_visits = True,

        # The epoch of coordinates for targets with measured proper motions. All catalogs must
        # match this epoch.
        epoch = 2016,

        # Generate full gurobi variable names instead of numbered ones (slow to build problem)
        # It is necessary only, when the netflow solution will be loaded from a file to
        # extract assignments based on variable names.
        use_named_variables = True,
    ),
    
    # These options are passed to the Gurobi LP solver and can significantly impact the time it
    # takes to find a solution
    gurobi_options = dict(
        # Random seed, set if you want the same results for the same input every time.
        # Do not set to avoid bad-case behavior of the solver.
        # seed = 0,

        # Agressiveness of presolver which tries to eliminate variables from the LP problem
        # before the simplex algorithm is run. -1 is automatic, 0 is off, 1-3 are different
        presolve = -1,

        # Method to use for solving the LP problem.
        # 3 means concurrent, 4 means deterministic concurrent
        method = 3,

        # Degenerate simplex moves, set to 0 to prevent too much time to be spent on
        # trying to improve the current solution
        degenmoves = 0,
        
        # How much of the time to spend by performing heuristics
        heuristics = 0.5,

        # mipfocus=1 is balanced toward finding more feasible solutions
        # mipfocus=2 is balanced toward proving that the current solution is the best
        # mipfocus=3 is to be used when the objection bound is moving very slowly
        mipfocus = 1,           
        
        # Relative stopping criterion for bounds on the objective. Set it to a few percent at least
        # because better solutions can take a very long time to find.
        mipgap = 0.01,

        # Generate log messages to the console
        LogToConsole = 1,
        
        # Time limit for the solver in seconds.
        timelimit = 300 # sec
    ),

    # The following debug options can be used to ignore certain constraints in the netflow
    # algorithm. These options are useful for debugging the algorithm, but should not be used
    # in production runs.
    debug_options = dict(
        # Ignore collisions between cobra tips
        ignore_endpoint_collisions = False,

        # Ignore collisions between cobra elbows
        ignore_elbow_collisions = False,

        # Ignore collisions between between broken cobras and their neighbors
        ignore_broken_cobra_collisions = False,

        # Ignore forbidden single targets
        ignore_forbidden_targets = False,

        # Ignore forbidden target pairs
        ignore_forbidden_pairs = False,

        # Ignore the minimum or the maximum number of fibers allocated to a target class
        ignore_calib_target_class_minimum = False,
        ignore_calib_target_class_maximum = False,

        # Ignore the minimum or the maximum number of fibers allocated to a target class
        ignore_science_target_class_minimum = False,
        ignore_science_target_class_maximum = False,

        # Ignore the time budget constraints
        ignore_time_budget = False,

        # Ignore the minimum or the maximum number of fibers in a cobra group allocated to a target class
        # If the number of calibration targets is low, ignoring the minimum requirement might be necessary
        # to produce a solution.
        ignore_cobra_group_minimum = False,
        ignore_cobra_group_maximum = False,

        # Ignores the `num_reserved_fibers` setting
        ignore_reserved_fibers = False,

        # Do not take the proper motion of stars into account when calculating the focal plane positions
        ignore_proper_motion = True,
    ),
)
