from collections.abc import Iterable
import numpy as np

def normalize_array(a, allow_none=True):
    if allow_none and a is None:
        return None
    elif a is None:
        raise ValueError('Argument cannot be None.')

    if isinstance(a, np.ndarray):
        return a

    if not isinstance(a, Iterable):
        a = [a]

    return np.array(a)