import os
import sys
import logging
import git
from git import GitCommandError
from datetime import datetime, timezone
from argparse import ArgumentParser
import commentjson as json
import yaml
import numpy as np

from ..util.notebookrunner import NotebookRunner

from ..setup_logger import logger

class Script():
    """
    Implements generic functions for command-line scripts.
    
    The class has support for attaching a debugger (debugpy + vscode) and
    profiling (cProfile). It also provides logging capabilities and argument
    parsing.
    """

    def __init__(self,
                 log_level=logging.INFO,
                 log_file=None,
                 log_to_file=True,
                 log_to_console=True):
        
        """
        Initialize the script with logging options.

        Parameters
        ---------
        log_level: int
            Logging level.
        log_file: str
            Log file name.
        log_to_file: bool
            Log to file.
        log_to_console: bool
            Log to console.
        """

        self.__log_level = log_level                # Default log level
        self.__log_to_file = log_to_file            # Log to file
        self.__log_to_console = log_to_console      # Log to console

        self.__debug = False                        # If True, script is running in debug mode
        self.__profile = False                      # If True, the profiler is enabled
        
        self.__log_formatter = None                 # Log formatter
        self.__log_file = log_file                  # Log file name
        self.__log_file_handler = None              # Log file handler
        self.__log_console_handler = None           # Log console handler

        self.__parser = ArgumentParser()
        self.__profiler = None
        self.__timestamp = datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')

    #region Properties

    def __get_debug(self):
        return self.__debug
    
    debug = property(__get_debug)

    def __get_profile(self):
        return self.__profile
    
    profile = property(__get_profile)

    def __get_log_level(self):
        return self.__log_level
    
    def __set_log_level(self, value):
        self.__log_level = value
    
    log_level = property(__get_log_level, __set_log_level)

    def __get_log_file(self):
        return self.__log_file
    
    def __set_log_file(self, value):
        self.__log_file = value
    
    log_file = property(__get_log_file, __set_log_file)

    def __get_log_to_file(self):
        return self.__log_to_file
    
    def __set_log_to_file(self, value):
        self.__log_to_file = value

    log_to_file = property(__get_log_to_file, __set_log_to_file)

    def __get_log_to_console(self):
        return self.__log_to_console
    
    def __set_log_to_console(self, value):
        self.__log_to_console = value

    log_to_console = property(__get_log_to_console, __set_log_to_console)

    def __get_timestamp(self):
        return self.__timestamp
    
    timestamp = property(__get_timestamp)

    #endregion

    def is_env(self, name):
        """
        Returns True if the environment variable `name` exists and it has a value.

        Parameters
        ---------
        name: str
            Environment variable name.

        Returns
        -------
        bool
            True if the environment variable exists and has a value.
        """

        return name in os.environ and os.environ[name] is not None
    
    def get_env(self, name, default=None):
        """
        Retuns the value of the environment variable `name` or `default` if it does not exist.

        Parameters
        ---------
            name: str
                Environment variable name.
            default: str
                Default value if the environment variable does not exist.

        Returns
        -------
        str
            Environment variable value or default.
        """

        if name in os.environ and os.environ[name] is not None:
            return os.environ[name]
        else:
            return default

    def __parse_args(self):
        """
        Parses the command-line arguments.
        """

        self.__args = self.__parser.parse_args().__dict__

    def add_arg(self, *args, **kwargs):
        """
        Adds an argument to the argument parser.

        Parameters
        ---------
        args: list
            Passed to `ArgumentParser.add_argument`.
        kwargs: dict
            Argument options.
        """

        self.__parser.add_argument(*args, **kwargs)

    def is_arg(self, name, args=None):
        """
        Returns True if the argument `name` exists and it has a value.

        Parameters
        ---------
        name: str
            Argument name.
        args: dict
            Arguments dictionary. If none, the arg dictionary from the last parse is used.

        Returns
        -------
        bool
            True if the argument exists and has a value.
        """

        args = args if args is not None else self.__args
        return name in args and args[name] is not None

    def get_arg(self, name, args=None, default=None):
        """
        Return the value of the argument `name` or `default` if it does not exist.

        Parameters
        ---------
        name: str
            Argument name.
        args: dict
            Arguments dictionary. If none, the arg dictionary from the last parse is used.
        default: any
            Default value if the argument does not exist.

        Returns
        -------
        any
            Argument value or default.
        """

        args = args if args is not None else self.__args

        if name in args and args[name] is not None and args[name] != '':
            return args[name]
        else:
            return default
        
    def _add_args(self):
        """
        Called by the script initialization routine to add command-line arguments
        from scripts. Override this method to add custom arguments.
        """

        self.add_arg('--debug', action='store_true', help='Enable debug mode')
        self.add_arg('--profile', action='store_true', help='Enable performance profiler')
        self.add_arg('--log-level', type=str, help='Set log level')

    def _init_from_args_pre_logging(self, args):
        """
        Initialize the most basic script settings from command-line arguments necessary
        to initialize the logging. Override this method to customize initialization.

        Parameters
        ---------
        args: dict
            Arguments dictionary.
        """

        self.__debug = self.get_arg('debug', args, self.__debug)
        self.__profile = self.get_arg('profile', args, self.__profile)

        if self.is_arg('log_level', args):
            log_level = self.get_arg('log_level', args)
            if isinstance(log_level, str) and hasattr(logging, log_level.upper()):
                self.__log_level = getattr(logging, log_level.upper())
            else:
                raise ValueError(f'Invalid log level `{log_level}`.')
            
        if self.__debug and self.__log_level > logging.DEBUG:
            self.__log_level = logging.DEBUG

    def _init_from_args(self, args):
        """
        Initialize script settings from command-line arguments. Override this method to
        customize initialization.

        Parameters
        ---------
        args: dict
            Arguments dictionary.
        """

        pass

    def _create_dir(self, name, dir, logger=logger):
        """
        Create a directory if it does not exist and log a message.

        Parameters
        ---------
        name: str
            Name of the directory.
        dir: str
            Directory path.
        logger: Logger
            Logger instance to be used. Default is the global logger.

        Returns
        -------
        bool
            True if the directory was created, False if it already exists.
        """

        dir = os.path.join(os.getcwd(), dir)
        if not os.path.isdir(dir):
            os.makedirs(dir, exist_ok=True)
            logger.debug(f'Created {name} directory `{dir}`.')
            return True
        else:
            logger.debug(f'Found existing {name} directory `{dir}`.')
            return False

    def get_command_name(self):
        """
        Returns the command name based on the script class name.

        Returns
        -------
        str
            Command name.
        """

        name = self.__class__.__name__.lower()
        name = name.replace('script', '')
        return f'ga-{name}'

    def start_logging(self):
        """
        Sets up logging for the script by configuring the root logger and adding
        handlers for console and file logging.
        """

        # Configure root logger
        root = logging.getLogger()
        root.handlers = []
        root.setLevel(self.__log_level)

        self.__log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s", datefmt='%H:%M:%S')
        
        if self.log_to_file and self.__log_file is not None:
            logdir = os.path.dirname(self.__log_file)
            self._create_dir('log', logdir)
            
            self.__log_file_handler = logging.FileHandler(self.__log_file)
            self.__log_file_handler.setFormatter(self.__log_formatter)

            root.addHandler(self.__log_file_handler)
        
        if self.__log_to_console:
            self.__log_console_handler = logging.StreamHandler()
            self.__log_console_handler.setFormatter(self.__log_formatter)

            root.addHandler(self.__log_console_handler)

        # Filter out log messages from matplotlib
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
 
        # Configure pipeline logger
        logger.propagate = True
        logger.setLevel(self.__log_level)

        if self.log_to_file and self.__log_file is not None:
            logger.info(f'Logging started to `{self.__log_file}`.')

    def suspend_logging(self):
        root = logging.getLogger()
        root.level = logging.CRITICAL

    def resume_logging(self):
        root = logging.getLogger()
        root.level = self.__log_level

    def stop_logging(self):
        """
        Stop logging and clean up log handlers.
        """

        if self.log_to_file and self.__log_file is not None:
            logger.info(f'Logging finished to `{self.__log_file}`.')

        # Disconnect file logger and re-assign stderr
        root = logging.getLogger()
        root.handlers = []
        root.addHandler(logging.StreamHandler())

        # Destroy logging objects (but keep last filename)
        self.__log_formatter = None
        self.__log_file_handler = None
        self.__log_console_handler = None

    def push_log_settings(self):
        """
        Stash the current log settings for later restoration.
        """

        root = logging.getLogger()
        self.__stashed_log_handlers = root.handlers

        self.__stashed_log_file_handler = self.__log_file_handler
        self.__stashed_log_formatter = self.__log_formatter
        self.__stashed_log_level = self.__log_level
        self.__stashed_log_file = self.__log_file

    def pop_log_settings(self):
        """
        Restore the stashed log settings.
        """

        root = logging.getLogger()
        root.handlers = self.__stashed_log_handlers

        self.__log_file_handler = self.__stashed_log_file_handler
        self.__log_formatter = self.__stashed_log_formatter
        self.__log_level = self.__stashed_log_level
        self.__log_file = self.__stashed_log_file

    def __start_profiler(self):
        """
        Start the profiler session
        """

        if self.__profile:
            import cProfile
            import pstats
            import io

            self.__profiler = cProfile.Profile()
            self.__profiler.enable()

            logger.info('Profiler started.')
        else:
            self.__profiler = None

    def __stop_profiler(self):
        """
        Stop the profiler session and save the results
        """

        if self.__profiler is not None:
            import pstats

            self.__profiler.disable()

            # Save profiler data to file
            with open(os.path.join('profile.cum.stats'), 'w') as f:
                ps = pstats.Stats(self.__profiler, stream=f).sort_stats('cumulative')
                ps.print_stats()

            with open(os.path.join('profile.tot.stats'), 'w') as f:
                ps = pstats.Stats(self.__profiler, stream=f).sort_stats('time')
                ps.print_stats()

            self.__profiler = None

            logger.info('Profiler stopped, results written to profile.*.stats.')

    def __dump_env(self, path):
        """
        Save environment variables to a file.        

        Parameters
        ---------
        path: str
            Full path to the file to save the environment variables.
        """

        with open(path, 'w') as f:
            for key, value in os.environ.items():
                f.write(f'{key}={value}\n')

        logger.debug(f'Environment variables saved to `{path}`.')

    def __dump_args(self, path, format=None):
        """
        Save the command-line arguments to a file. The arguments that are
        saved are the effective ones, ie. they might have been modified
        by the script.

        Parameters
        ---------
        path: str
            Full path to the file to save the arguments.
        """

        def default(obj):
            if isinstance(obj, float):
                return "%.5f" % obj
            if type(obj).__module__ == np.__name__:
                if isinstance(obj, np.ndarray):
                    if obj.size < 100:
                        return obj.tolist()
                    else:
                        return "(not serialized)"
                else:
                    return obj.item()
            return "(not serialized)"

        if format is None:
            _, ext = os.path.splitext(path)
            with open(path, 'w') as f:
                if ext in [ '.json', '.yaml' ]:
                    format = ext

        with open(path, 'w') as f:
            if format == '.json':
                json.dump(self.__args, f, default=default, indent=4)
            elif format == '.yaml':
                yaml.dump(self.__args, f, indent=4)

        logger.debug(f'Arguments saved to `{path}`.')
            
    def __dump_cmdline(self, path):
        """
        Save to command-line into a file.

        Parameters
        ---------
        path: str
            Full path to the file to save the command-line.
        """

        with open(path, 'w') as f:
            f.write(f'{self.get_command_name()} ')
            if len(sys.argv) > 1:
                f.write(' '.join(sys.argv[1:]))
            f.write('\n')

        logger.debug(f'Command-line saved to `{path}`.')

    def _dump_settings(self):
        """
        Save environment, arguments and command-line to files next to the log file
        """

        if self.__log_to_file and self.__log_file is not None:
            logdir = os.path.dirname(self.__log_file)
            command = self.get_command_name()
            self.__dump_env(os.path.join(logdir, f'{command}_{self.__timestamp}.env'))
            self.__dump_args(os.path.join(logdir, f'{command}_{self.__timestamp}.args'), format='.json')
            self.__dump_cmdline(os.path.join(logdir, f'{command}_{self.__timestamp}.cmd'))

    def execute(self):
        """
        Execute the script by calling the life-cycle functions in order.

        Do not override this method.
        """

        # Initialize the command-line parser and parse the arguments
        self._add_args()
        self.__parse_args()

        # Initialize the script settings before logging is started
        self._init_from_args_pre_logging(self.__args)

        # NOTE: debugging is started from the wrapper shell script

        # Initialize output directory, override log file location, etc.
        self.prepare()

        # Start logging and profiler
        self.start_logging()    
        self._init_from_args(self.__args)
        self._dump_settings()
        self.__start_profiler()

        logger.info(f'Starting execution of {self.get_command_name()}.')

        try:
            self.run()
        except Exception as ex:
            # Log the exception first, it will write the stack trace to the log
            logger.exception('Unhandled exception occured.')

            # If running in the debugger, break here so that we can go back and
            # check the error
            try:
                import debugpy
                if debugpy.is_client_connected:
                    raise ex
            except ImportError:
                pass

        logger.info(f'Finished execution of {self.get_command_name()}.')

        self.__stop_profiler()
        self.stop_logging()

    def prepare(self):
        """
        Executed before the main run method. 
        
        Override this method and do initializations here, such as setting
        up the logging level, directories, etc.
        """
        
        command = self.get_command_name()
        time = self.__timestamp
        self.__log_file = f'{command}_{time}.log'

    def run(self):
        """
        Override this function to implement the main logic of the script.
        """
        
        raise NotImplementedError()
    

    def get_last_git_commit(self, module):
        dir = os.path.dirname(os.path.abspath(module.__file__))
        repo = git.Repo(dir, search_parent_directories=True)

        current_hash = repo.head.object.hexsha

        try:
            recent_tag = repo.git.describe(tags=True, abbrev=0)
        except GitCommandError:
            recent_tag = None

        try:
            current_branch = repo.active_branch.name
        except TypeError:
            current_branch = None

        return current_hash, current_branch, recent_tag

    def _execute_notebook(self, notebook_path, parameters, outdir):
        fn = os.path.basename(notebook_path)
        logger.info(f'Executing notebook `{fn}`.')

        nr = NotebookRunner()
        nr.workdir = os.path.dirname(notebook_path)
        nr.parameters = parameters
        # nr.kernel = kernel
        nr.open_ipynb(notebook_path)

        # Suspend logging because nbconvert writes and insane
        # amount of messages at debug level
        self.suspend_logging()
        nr.run()
        self.resume_logging()

        fn, ext = os.path.splitext(os.path.basename(notebook_path))
        nr.save_ipynb(os.path.join(outdir, fn + '.ipynb'))
        nr.save_html(os.path.join(outdir, fn + '.html'))