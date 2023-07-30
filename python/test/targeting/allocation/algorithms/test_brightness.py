import os
import numpy as np
import matplotlib.pyplot as plt

from ics.cobraOps.Bench import Bench
from pfs.ga.targeting.allocation import algorithms

from test_base import TestBase

from pfs.ga.targeting.projection import Pointing
from pfs.ga.targeting.instrument import SubaruPFI, SubaruWFC
from pfs.ga.targeting.allocation.algorithms import Brightness

class BrightnessTest(TestBase):
    def test_argsort_assoc(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()   
        wfc = SubaruWFC(pointing=Pointing(227.2420, 67.2221))
        pfi = SubaruPFI(projection=wfc)

        cmd, _ = self.get_test_cmd()
        algorithm = Brightness(cmd.axes[1].magnitude)

        mask = np.full_like(ra, True, dtype=bool)
        mask[1000:] = False

        assoc = pfi.find_associations(obs, mask=mask)
        sorting = algorithm.argsort_assoc(obs, assoc, mask=mask)

        ix, ixm = assoc.pick_first(sorting)

        pass

    def test_combine_assoc(self):
        obs = self.load_test_observation()
        ra, dec = obs.get_coords()

        mask = np.full(ra.shape, True)
        targets = np.full(ra.shape, False)

        for pp in [(228.2, 67.5, 0),
                   (226.3, 67.5, 0),
                   (226.0, 66.9, 0),
                   (228.1, 66.955, 40)]:

            wfc = SubaruWFC(pointing=Pointing(*pp))
            pfi = SubaruPFI(projection=wfc)

            cmd, _ = self.get_test_cmd()
            algorithm = Brightness(cmd.axes[1].magnitude)

            assoc = pfi.find_associations(obs, mask=mask)
            sorting = algorithm.argsort_assoc(obs, assoc, mask=mask)
            ix, ixm = assoc.pick_first(sorting)

            midx = tuple(a[ix[ixm]] for a in np.where(mask))

            mask[midx] = False
            targets[midx] = True

        pass
