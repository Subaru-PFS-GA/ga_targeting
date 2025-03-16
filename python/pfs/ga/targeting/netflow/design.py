import os
import numpy as np
from collections.abc import Iterable

from pfs.ga.targeting.instrument import SubaruPFI

class Design():
    """
    Utility class to create PfsDesign object from a list of targets.
    """

    def __init__(self):
         pass

    def join_catalogs(assignments, catalogs):
        catalogs = [catalogs] if not isinstance(catalogs, Iterable) else catalogs

        for catalog in catalogs:
            assignments = assignments.merge(catalog, on='targetid', how='left')
    