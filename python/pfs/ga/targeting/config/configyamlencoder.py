import yaml
import numpy as np
from datetime import datetime, date

def config_yaml_array_representer(dumper, data):
    return dumper.represent_list(data.tolist())

def config_yaml_scalar_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:float', data.item())

def config_yaml_datetime_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:timestamp', data.isoformat())

def config_yaml_date_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:timestamp', data.isoformat())

yaml.add_representer(np.ndarray, config_yaml_array_representer)
yaml.add_representer(np.generic, config_yaml_scalar_representer)
yaml.add_representer(datetime, config_yaml_datetime_representer)
yaml.add_representer(date, config_yaml_date_representer)