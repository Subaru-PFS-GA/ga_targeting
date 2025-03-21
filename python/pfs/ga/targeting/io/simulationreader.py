from ..data import Simulation
from .catalogserializer import CatalogSerializer

class SimulationReader(CatalogSerializer):
    def __init__(self, orig=None):
        super(SimulationReader, self).__init__(orig=orig)

        if not isinstance(orig, SimulationReader):
            pass
        else:
            pass

    def _create_catalog(self, name=None):
        sim = Simulation(name=name)
        for _, p in self.photometry.items():
            sim.append_photometry(p)
        return sim