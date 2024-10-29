import os
import sys
import logging
from datetime import datetime, timezone
from argparse import ArgumentParser
import commentjson as json
import yaml
import numpy as np

from ..setup_logger import logger

class Script():
    """
    Implements generic function for pipeline command-line scripts.
    """

    def __init__(self,
                 log_level=logging.INFO,
                 log_file=None,
                 log_to_file=True,
                 log_to_console=True):

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
        return name in os.environ and os.environ[name] is not None
    
    def get_env(self, name, default=None):
        if name in os.environ and os.environ[name] is not None:
            return os.environ[name]
        else:
            return default

    def __parse_args(self):
        self.__args = self.__parser.parse_args().__dict__

    def add_arg(self, *args, **kwargs):
        self.__parser.add_argument(*args, **kwargs)

    def is_arg(self, name, args=None):
        args = args if args is not None else self.__args
        return name in args and args[name] is not None

    def get_arg(self, name, args=None, default=None):
        args = args if args is not None else self.__args

        if name in args and args[name] is not None and args[name] != '':
            return args[name]
        else:
            return default
        
    def _add_args(self):
        self.add_arg('--debug', action='store_true', help='Enable debug mode')
        self.add_arg('--profile', action='store_true', help='Enable performance profiler')
        self.add_arg('--log-level', type=str, help='Set log level')

    def _init_from_args(self, args):
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

    def _create_dir(self, name, dir, logger=logger):
        dir = os.path.join(os.getcwd(), dir)
        if not os.path.isdir(dir):
            os.makedirs(dir, exist_ok=True)
            logger.debug(f'Created {name} directory `{dir}`.')
        else:
            logger.debug(f'Found existing {name} directory `{dir}`.')

    def get_command_name(self):
        return f'targeting-{self.__class__.__name__.lower()}'

    def start_logging(self):
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

    def stop_logging(self):
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

        root = logging.getLogger()
        root.handlers = self.__stashed_log_handlers

        self.__log_file_handler = self.__stashed_log_file_handler
        self.__log_formatter = self.__stashed_log_formatter
        self.__log_level = self.__stashed_log_level
        self.__log_file = self.__stashed_log_file

    def pop_log_settings(self):
        """
        Restore the stashed log settings.
        """

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

    def __dump_env(self, filename):
        with open(filename, 'w') as f:
            for key, value in os.environ.items():
                f.write(f'{key}={value}\n')

        logger.debug(f'Environment variables saved to `{filename}`.')

    def __dump_args(self, filename):

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

        _, ext = os.path.splitext(filename)
        with open(filename, 'w') as f:
            if ext == '.json':
                json.dump(self.__args, f, default=default, indent=4)
            elif ext == '.yaml':
                yaml.dump(self.__args, f, indent=4)

        logger.debug(f'Arguments saved to `{filename}`.')
            
    def __dump_cmdline(self, filename):
        """
        Save to command-line into a file.
        """

        with open(filename, 'w') as f:
            f.write(f'{self.get_command_name()} ')
            if len(sys.argv) > 1:
                f.write(' '.join(sys.argv[1:]))
            f.write('\n')

        logger.debug(f'Command-line saved to `{filename}`.')

    def _dump_settings(self):
        # Save environment, arguments and command-line to files next to the log file
        if self.__log_to_file and self.__log_file is not None:
            logdir = os.path.dirname(self.__log_file)
            command = self.get_command_name()
            self.__dump_env(os.path.join(logdir, f'{command}_{self.__timestamp}.env'))
            self.__dump_args(os.path.join(logdir, f'{command}_{self.__timestamp}.args.json'))
            self.__dump_cmdline(os.path.join(logdir, f'{command}_{self.__timestamp}.cmd'))

    def execute(self):
        self._add_args()
        self.__parse_args()
        self._init_from_args(self.__args)

        # NOTE: debugging is started from the wrapper shell script

        self.prepare()

        self.start_logging()    
        self._dump_settings()
        self.__start_profiler()

        self.run()

        self.__stop_profiler()
        self.stop_logging()

    def prepare(self):
        """
        Executed before the main run method. Do initializations here, such as setting
        up the logging level, directories, etc.
        """
        
        command = self.get_command_name()
        time = self.__timestamp
        self.__log_file = f'{command}_{time}.log'

    def run(self):
        raise NotImplementedError()