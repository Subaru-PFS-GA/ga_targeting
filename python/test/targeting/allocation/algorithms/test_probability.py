import os
import numpy as np
import matplotlib.pyplot as plt

from ics.cobraOps.Bench import Bench
from pfs.ga.targeting.allocation import algorithms
from pfs.ga.targeting.allocation.algorithms.probability import Probability

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC
from pfs.ga.targeting.allocation.algorithms import Brightness
from pfs.ga.targeting import ProbabilityMap

class ProbabilityTest(TestBase):
    def test_sort_assoc(self):
        obs = self.load_test_observation()
        sim = self.load_test_simulation()
        cmd, photometry = self.get_test_cmd()
        
        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)
        pmap.merge_populations((np.s_[:10], np.s_[10:]))

        lp_member, mask_member = pmap.lookup_lp_member(obs)
        
        wfc = SubaruWFC(pointing=Pointing(227.2420, 67.2221))        
        pfi = SubaruPFI(projection=wfc)
        assoc = pfi.find_associations(obs)

        algorithm = Probability()
        ix = algorithm.argsort_assoc(lp_member, assoc)
        ix, mask = assoc.pick_first(ix)

        pass