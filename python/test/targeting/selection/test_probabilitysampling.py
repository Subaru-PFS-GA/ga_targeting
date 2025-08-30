import numpy as np

from pfs.ga.common.photometry import Photometry, Magnitude, Color
from pfs.ga.targeting import ProbabilityMap
from pfs.ga.targeting.selection import ProbabilitySampling

from test_base import TestBase

class ProbabilitySamplingTest(TestBase):
    def test_apply(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()
        obs = self.load_test_observation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)
        pmap.merge_populations((np.s_[:10], np.s_[10:]))

        sel = ProbabilitySampling(pmap, 1)
        mask = sel.apply(obs)

        pass