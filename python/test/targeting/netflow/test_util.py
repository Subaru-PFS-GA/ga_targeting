import os

from pfs.ga.common.util.config import *
from test_base import TestBase

class InstrumentTest(TestBase):
    def test_camel_to_snake(self):
        self.assertIsNone(camel_to_snake(None))
        self.assertEqual(camel_to_snake('Camel'), 'camel')
        self.assertEqual(camel_to_snake('CamelCase'), 'camel_case')
        self.assertEqual(camel_to_snake('CamelCaseCase'), 'camel_case_case')
        self.assertEqual(camel_to_snake('camel'), 'camel')
        self.assertEqual(camel_to_snake('camelCase'), 'camel_case')
        self.assertEqual(camel_to_snake('camelCaseCase'), 'camel_case_case')