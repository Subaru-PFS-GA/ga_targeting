from pfs.ga.targeting.scripts.targetingscript import TargetingScript

from test_base import TestBase

class TargetingScriptTest(TestBase):
    def test_init(self):
        script = TargetingScript()
        self.assertIsNotNone(script)
