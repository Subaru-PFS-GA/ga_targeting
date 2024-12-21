import os
import numpy as np
import commentjson as json
import yaml
from typing import List, Dict, get_type_hints, get_origin, get_args

from ..setup_logger import logger
from .configjsonencoder import ConfigJSONEncoder
from .configyamlencoder import *

class Config():
    """
    Base class for configurations objects. Implements function to convert
    a hierarchy of configuration classes into a dictionary and back. The
    dictionaries, in turn, can be saved to a file or loaded from a file. It
    supports loading dictionaries from Python files (via `exec`), JSON files
    and YAML files.

    The class can keep track of the configuration files used to load the configuration
    from.
    """

    def __init__(self):
        self.__config_files = None                # List of configuration files loaded

    @classmethod
    def default(cls):
        """
        Create a default configuration.
        """

        return cls()
    
    def to_dict(self):
        """
        Convert the configuration to a dictionary.
        """

        return self._save_impl()
    
    def __repr__(self):
        return str(self.to_dict())

    #region Properties

    def __get_config_files(self):
        return self.__config_files

    config_files = property(__get_config_files)

    #endregion
    #region Utility functions

    def _get_env(self, name, default=None):
        if name in os.environ and os.environ[name] is not None and os.environ[name] != '':
            return os.environ[name]
        else:
            return None
        
    #endregion
    #region Load

    @classmethod
    def from_file(cls, path, format=None, ignore_collisions=False):
        """
        Load the configuration from a file.

        Arguments
        ---------
        path : str
            Path to the configuration file.
        ignore_collisions : bool
            If True, collisions in the configuration are ignored. If False, an exception is
            raised when a collision is detected.
        """

        c = cls()
        c.load(path, format=format, ignore_collisions=ignore_collisions)
        return c
    
    @classmethod
    def from_dict(cls, config, ignore_collisions=False):
        """
        Load the configuration from a dictionary.

        Arguments
        ---------
        config : dict
            Dictionary containing the configuration.
        ignore_collisions : bool
            If True, collisions in the configuration are ignored. If False, an exception is
            raised when a collision is detected.
        """

        c = cls()
        c.load(config, ignore_collisions=ignore_collisions)
        return c

    def load(self, source, format=None, ignore_collisions=False):
        """
        Load the configuration from a dictionary or file.

        Arguments
        ---------
        source : dict, str, list
            Source of the configuration. If a dictionary, the configuration is loaded from
            the dictionary. If a string, the configuration is loaded from a file. If a list
            of strings, the configurations are loaded from the list of files.
        ignore_collisions : bool
            If True, collisions in the configuration are ignored. If False, an exception is
            raised when a collision is detected.
        """

        if source is None:
            pass
        elif isinstance(source, dict):
            # Configuration is a dictionary, no need to read from a file
            self._load_impl(config=source, ignore_collisions=ignore_collisions)
        elif isinstance(source, str) or isinstance(source, list):
            # Configuration to be loaded from a single file or from a list of files
            # When loading from a list of files, the configuration is merged

            # Keep track of config files used
            if self.__config_files is None:
                self.__config_files = []

            source = source if isinstance(source, list) else [ source ]
            for s in source:
                # Load the source file as a dictionary
                config = Config.__load_dict_from_file(s, format=format)
                
                # Load the configuration from the dictionary
                self._load_impl(config=config, ignore_collisions=ignore_collisions)

                # Keep track of the configuration files used
                self.__config_files.append(os.path.abspath(s))
        else:
            raise NotImplementedError()
    
    def _load_impl(self, config=None, ignore_collisions=False):
        """
        Config class specific implementation of how to load the configuration from
        a dictionary. Override this function in subclasses to implement custom
        functionality.
        """

        # The default is not to use type mapping
        if config is not None:
            self.__load_config_from_dict(config=config, ignore_collisions=ignore_collisions)
    
    def __load_config_from_dict(self, config=None, type_hints=None, ignore_collisions=False):
        """
        Load configuration from a dictionary.

        This function iterates over the dictionary passed as `config` and sets the members
        of the configuration class to the values found in the dictionary. If a member is a subclass
        of `Config`, the value found in the dictionary (usually a sub-dictionary) is passed to the class
        for further processing. If key in the dictionary matches a key in `type_map`, then the value
        is passed to the class defined in the map.
        
        If the member is a dictionary and the value in the dictionary is also a dictionary, the two
        dictionaries are merged. If the member is a dictionary and the value in the dictionary is not
        a dictionary, the value is set as is. If the member is not a dictionary, the value is set as is.
        """

        annotations = type(self).__init__.__annotations__

        # Iterate over all keys of the configuration and see if the keys match up with
        # member variables of the configuration class
        for key, value in config.items():
            if not hasattr(self, key):
                raise ValueError(f'Member `{key}` of class `{type(self).__name__}` does not exist.')
            
            if value is None:
                setattr(self, key, None)
                continue
            
            # If the member is found and it's a subclass of `Config`, just pass the dictionary
            # to it for further processing. If the member is found but its value if not a subclass
            # of `Config` but its name is in `type_map`, then instantiate the particular type defined
            # in the map. In all other cases, just set the member to the value found in the config dict.

            c = getattr(self, key)

            if isinstance(c, Config):
                # This is a config class, it can initialize itself from the config dict
                c._load_impl(value, ignore_collisions=ignore_collisions)
            elif key in annotations:
                if isinstance(annotations[key], type) and issubclass(annotations[key], Config):
                    # This member is a type-annotated member where the type is a subclass of Config,
                    # convert the value to the type
                    setattr(self, key, Config.__config_to_class(annotations[key], config=value, existing=None, ignore_collisions=ignore_collisions))
                else:
                    # This is an annotation in the form of List or Dict or similar
                    args = get_args(annotations[key])
                    if any([ issubclass(a, Config) for a in args ]):
                        # This is a list or dictionary of Config classes
                        setattr(self, key, Config.__config_to_class(annotations[key], config=value, existing=c, ignore_collisions=ignore_collisions))
                    else:
                        # This is a regular member with some annotations, just set its value
                        setattr(self, key, value)
            elif type_hints is not None and key in type_hints:
                # This member is part of the type_map, instantiate the type and initialize
                setattr(self, key, Config.__config_to_class(type_hints[key], config=value, existing=c, ignore_collisions=ignore_collisions))
            elif isinstance(c, dict) and isinstance(value, dict):
                # This member is a dictionary, merge with the config dict
                setattr(self, key, self.__merge_dict(c, value, ignore_collisions=ignore_collisions))
            else:
                # This is a regular member, just set its value
                setattr(self, key, value)

    @staticmethod
    def __config_to_class(class_type, config=None, existing=None, ignore_collisions=False):
        """
        Convert a configuration to a class instance. The type of the target class
        is passed as type hint which can also be a list or dict.

        Parameters
        ----------
        class_type : type
            Type hint of the target class.
        config : dict
            Configuration to be converted to a class instance.
        existing : object
            Existing object to be updated with the configuration.
        """

        # Check if the type hint is a list or a dictionary and get the type of the elements
        origin = get_origin(class_type)
        args = get_args(class_type)
        
        if origin == list:
            # This is a list of config entries that need to be converted to well-known config classes
            # We can't verify which object is which, so the best we can do is to fully replace the list
            return [ Config.__config_to_class(args[0], c, ignore_collisions=ignore_collisions) for c in config ]
        elif origin == dict:
            # This is a dictionary to be mapped to well-known config classes
            # Iterate through the dictionary and update each class, if already there,
            # otherwise create a new instance
            # TODO: how do we force creating new and removing items?
            if existing is None:
                existing = {}

            for k, c in config.items():
                if k in existing:
                    # Update existing config object
                    existing[k]._load_impl(config=c, ignore_collisions=ignore_collisions)
                else:
                    # Create a new config object
                    existing[k] = Config.__config_to_class(args[1], c, ignore_collisions=ignore_collisions)

            return existing
        else:
            v = class_type()
            v.load(config, ignore_collisions=ignore_collisions)        
            return v

    @staticmethod
    def __load_dict_from_file(path, format=None):
        """
        Depending on the file extension, load the configuration file
        """

        if format is None:
            dir, filename = os.path.split(path)
            _, format = os.path.splitext(filename)

        if format == '.py':
            config = Config.__load_dict_py(path)
        elif format == '.json':
            config = Config.__load_dict_json(path)
        elif format == '.yaml':
            config = Config.__load_dict_yaml(path)
        else:
            raise ValueError(f'Unknown configuration file extension `{format}`')
        
        return config

    @staticmethod
    def __load_dict_py(filename):
        # Load a python file with the configuration and execute it to get
        # the configuration dictionary

        with open(filename, 'r') as f:
            code = f.read()

        # TODO: something is wrong here because lambda functions using
        #       numpy wont work, maybe we should pass global() and local()?

        global_variables = {}
        local_variables = {}
        exec(code, global_variables, local_variables)

        if 'config' in local_variables:
            return local_variables['config']
        else:
            raise ValueError(f'Configuration not found in file `{filename}`')
    
    @staticmethod
    def __load_dict_json(filename):
        # Load configuration from a JSON file with comments

        with open(filename, 'r') as f:
            config = json.load(f)
        return config
    
    @staticmethod
    def __load_dict_yaml(filename):
        # Load configuration from a YAML file

        with open(filename, 'r') as f:
            config = yaml.safe_load(f)
        return config

    #endregion Load
        
    def save(self, path):
        # Save configuration to a file

        config = self._save_impl()
        Config.__save_dict_to_file(config, path, format='.json')

    def _save_impl(self):
        """
        Config class specific implementation of how to save the configuration to
        a dictionary. Override this function in subclasses to implement custom
        functionality.
        """

        return self.__save_config_to_dict(self)
    
    @staticmethod
    def __save_config_to_dict(obj):
        """
        Save configuration to a dictionary
        """

        config = {}
        for k in obj.__dict__:
            # Save public members only
            if not k.startswith('_'):
                v = getattr(obj, k)
                config[k] = Config.__class_to_config(v)
        return config
    
    @staticmethod
    def __class_to_config(obj):
        """
        Convert a configuration class instance to a dictionary by iterating over
        its public members. In addition, list and dictionaries are kept but
        numpy arrays are unwrapped into basic python lists and numpy numbers are
        casted to the closes python types.
        """

        if isinstance(obj, Config):
            return Config.__save_config_to_dict(obj)
        elif isinstance(obj, np.ndarray):
            return [ Config.__class_to_config(v) for v in obj.tolist() ]
        elif isinstance(obj, np.generic):
            return obj.item()
        elif isinstance(obj, dict):
            return { k: Config.__class_to_config(v) for k, v in obj.items() }
        elif isinstance(obj, list):
            return [ Config.__class_to_config(v) for v in obj ]
        else:
            return obj
          
    @staticmethod
    def __save_dict_to_file(config, path, format=None):
        """
        Depending on the file extension, save the configuration file
        """

        dir, filename = os.path.split(path)

        if format is None:
            _, ext = os.path.splitext(filename)
            if ext in [ '.py', '.json', '.yaml' ]:
                format = ext
            else:
                raise ValueError(f'Unknown configuration file extension `{ext}`')
        
        if format == '.py':
            Config.__save_dict_py(config, path)
        if format == '.json':
            Config.__save_dict_json(config, path)
        elif format == '.yaml':
            Config.__save_dict_yaml(config, path)
        else:
            raise ValueError(f'Unknown file format `{format}`')

    @staticmethod
    def __save_dict_py(config, filename):
        # Save a python file with the configuration
        raise NotImplementedError()

    @staticmethod
    def __save_dict_json(config, filename):
        # Save configuration to a JSON file with comments
        with open(filename, 'w') as f:
            json.dump(config, f,
                      sort_keys=False,
                      indent=2,
                      cls=ConfigJSONEncoder)

    @staticmethod
    def __save_dict_yaml(config, filename):
        # Save configuration to a YAML file
        with open(filename, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    #region Dictionary utilities
        
    @staticmethod
    def __merge_dict(a: dict, b: dict, ignore_collisions=False):
        """
        Deep-merge two dictionaries. This function will merge the two dictionaries
        recursively. If a key is present in both dictionaries, the value will be 
        merged recursively, unless a collision is detected.

        Since file formats such as yaml and json do not support all data types,
        if dictionary `a` already contains a key, all keys in `b` are converted to
        the same data type.
        """

        # Verify that all keys in dictionary `a` have the same type
        type_a = None
        for k in a.keys():
            if type_a is None:
                type_a = type(k)
            elif type(k) != type_a:
                raise ValueError(f"Key `{k}` in dictionary `a` has a different type than the other keys.")

        # If dictionary `a` has a key type, use that, otherwise try to determine the key type for
        # dictionary `b` starting with int, float and finally falling back to str
        bb = None
        if type_a is not None:
            try:
                bb = { type_a(k): b[k] for k in b.keys() }
            except ValueError:
                raise ValueError(f"Key type mismatch between dictionary `a` and `b`.")
        else:
            for t in [ int, float, str ]:
                try:
                    bb = { t(k): b[k] for k in b.keys() }
                    break
                except ValueError:
                    pass

        kk = list(a.keys()) + list(bb.keys())

        r = {}
        for k in kk:
            # Both are dictionaries, merge them
            if k in a and isinstance(a[k], dict) and k in bb  and isinstance(bb[k], dict):
                r[k] = Config.__merge_dict(a[k], bb[k], ignore_collisions=ignore_collisions)
            elif k in a and k in bb:
                msg = f"Collision detected in the configuration for key `{k}`."
                if ignore_collisions:
                    r[k] = bb[k]
                    logger.debug(msg)
                else:
                    raise ValueError(msg)
            elif k in a:
                r[k] = a[k]
            elif k in bb:
                r[k] = bb[k]

        return r
    
    @staticmethod
    def __copy_dict(a: dict):
        # Make a deep copy of a dictionary
        r = {}
        for k in a.keys():
            if isinstance(a[k], dict):
                r[k] = Config.__copy_dict(a[k])
            elif isinstance(a[k], list):
                r[k] = [ Config.__copy_dict(c) for c in a[k] ]
            else:
                r[k] = a[k]
        return r
    
    #endregion

