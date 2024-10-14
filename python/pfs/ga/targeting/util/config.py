import re

def camel_to_snake(s):
    """
    If a dictionary (of settings) is provided with camelCase keys, convert them to snake_case.
    """

    if s is None:
        return None
    elif isinstance(s, dict):
        return { camel_to_snake(k): v for k, v in s.items() }
    elif isinstance(s, list):
        return s
    elif isinstance(s, str):
        if '_' in s:
            return s
        else:
            return re.sub(r'(?<!^)(?=[A-Z])', '_', s).lower()
    else:
        raise NotImplementedError()