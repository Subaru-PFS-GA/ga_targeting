import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

import pfs.ga.targeting
from pfs.ga.targeting.targets.dsph import GALAXIES as DSPH_FIELDS

class NetflowConfigTest(TestBase):
    def test_get_netflow_config(self):
        field = DSPH_FIELDS['sextans']
        config = field.get_netflow_config()
