import os
import numpy as np

from pfs.ga.common.photometry import Photometry, Magnitude, Color
from pfs.ga.common.selection import MagnitudeSelection
from pfs.ga.targeting import ProbabilityMap

from test_base import TestBase

class ProbabilityMapTest(TestBase):
    def test_from_simulation(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        # Without mask
        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)

        # With magnitude mask
        mask = MagnitudeSelection(cmd.axes[1], 17, 22.5).apply(sim)
        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim, mask=mask)

    def test_merge_populations(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)

        pmap.merge_populations((np.s_[:10], np.s_[10:]))

    def test_maximum_filter(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)
        
        pmap.maximum_filter()

    def test_get_nonzero_mask(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)

        pmap.maximum_filter()
        pmap.merge_populations((np.s_[:10], np.s_[10:]))

        mask = pmap.get_nonzero_mask()
        mask = pmap.get_nonzero_mask(populations=[0])

    def test_get_norm(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)

        pmap.maximum_filter()
        pmap.merge_populations((np.s_[:10], np.s_[10:]))
        mask = pmap.get_nonzero_mask()

        norm = pmap.get_norm()
        self.assertEqual((2,), norm.shape)

        norm = pmap.get_norm(mask=mask)
        self.assertEqual((2,), norm.shape)

    def test_get_lp_member(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim, bins=[100, 100], merge_list=[np.s_[:10], np.s_[10:]])
        pmap.maximum_filter()
        
        
        lp_member, mask_member = pmap.get_lp_member()
        self.assertEqual((2, 100, 100), lp_member.shape)

    def test_lookup_lp_member(self):
        cmd, photometry = self.get_test_cmd()
        sim = self.load_test_simulation()
        obs = self.load_test_observation()

        mask = MagnitudeSelection(cmd.axes[1], 17, 22.5).apply(obs)

        pmap = ProbabilityMap(cmd.axes)
        pmap.from_simulation(sim)
        pmap.merge_populations((np.s_[:10], np.s_[10:]))

        # Without per-star weights
        lp_member, lp_member_mask = pmap.lookup_lp_member(obs)
        lp_member, lp_member_mask = pmap.lookup_lp_member(obs, mask=mask)

        # With per-star weights
        w = np.random.dirichlet([0.5, 0.5], size=obs.shape)
        lp_member, lp_member_mask = pmap.lookup_lp_member(obs, w)
        lp_member, lp_member_mask = pmap.lookup_lp_member(obs, w[mask], mask=mask)