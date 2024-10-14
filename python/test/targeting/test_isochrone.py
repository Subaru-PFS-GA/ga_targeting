import os

from test_base import TestBase
from pfs.ga.targeting.photometry import Photometry, Magnitude
from pfs.ga.targeting.diagram import CMD
from pfs.ga.targeting.instrument import SubaruHSC
from pfs.ga.targeting import Isochrone

class IsochroneTest(TestBase):
    def test_from_isogrid(self):
        iso = self.get_test_isochrone()

    def test_get_magnitude(self):
        iso = self.get_test_isochrone()    
        cmd, photometry = self.get_test_cmd()
        m1, s1 = iso.get_magnitude(cmd.axes[1].magnitude, DM=19)
        
    def test_get_color(self):
        iso = self.get_test_isochrone()
        cmd, photometry = self.get_test_cmd()
        c1, s1 = iso.get_color(cmd.axes[0].color)
