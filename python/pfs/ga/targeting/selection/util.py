from collections.abc import Iterable
import numpy as np

def validate_photometry(axes):
    """
    Make sure axes use the same photometric system.
    """

    p = None
    for a in axes:
        if p is None:
            p = a.photometry
        elif p != a.photometry:
            raise Exception('Photometric systems must match.')