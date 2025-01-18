import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

import pfs.ga.targeting
from pfs.ga.targeting.config import NetflowConfig

class NetflowConfigTest(TestBase):
    def test_init(self):
        config = NetflowConfig()

    def test_default(self):
        config = NetflowConfig.default()

    def test_example(self):
        path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../configs/netflow/example.py')
        config = NetflowConfig.from_file(path)

        self.assertEqual(config.field.nvisits, 6)
        
        self.assertEqual(len(config.pointings), 2)
        
        self.assertEqual(len(config.targets), 3)
        self.assertEqual(len(config.targets['dsph'].filters), 2)
        self.assertIsNone(config.targets['dsph'].bands)
        self.assertIsNone(config.targets['fluxstd'].filters)
        self.assertEqual(len(config.targets['fluxstd'].bands), 2)
        
        self.assertIsNotNone(config.instrument_options)

        dist = np.random.rand(10)
        penalty = config.netflow_options.black_dot_penalty(dist)

        self.assertEqual(len(config.netflow_options.target_classes), 8)
        self.assertEqual(len(config.netflow_options.cobra_groups), 2)
        
        self.assertEqual(len(config.netflow_options.time_budgets), 1)
        self.assertEqual(len(config.netflow_options.time_budgets['science'].target_classes), 6)

        self.assertIsNotNone(config.gurobi_options)

        self.assertIsNotNone(config.debug_options)

        path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../tmp/netflow/example.json')
        config.save(path)