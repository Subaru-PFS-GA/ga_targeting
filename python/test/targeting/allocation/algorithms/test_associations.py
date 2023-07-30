import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC

class AssociationsTest(TestBase):
    def test_find_multiplets(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()
        mask = np.full_like(ra, True, dtype=bool)
        
        wfc = SubaruWFC(pointing=Pointing(227.2420, 67.2221))
        pfi = SubaruPFI(projection=wfc)
        assoc = pfi.find_associations(obs, mask=mask)

        ix, counts = assoc.find_multiplets()
        pass
        