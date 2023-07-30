import h5py

from .simulationreader import SimulationReader

class Hdf5SimulationReader(SimulationReader):
    def __init__(self, orig=None):
        super(Hdf5SimulationReader, self).__init__(orig=orig)

        if not isinstance(orig, Hdf5SimulationReader):
            pass
        else:
            pass

    def read(self, filename, mask=None, name=None, **kwargs):

        # TODO: this is a bit simplified, we could read all physical parameters etc.
        # d = {}
        # with h5py.File(filename, 'r') as h:
        #     for _, p in self.photometry.items():
        #         for _, m in p.magnitudes.items():
        #             for prefix in ['', 'obs_', 'err_']:
        #                 k = m.get_name(prefix=prefix)
        #                 if k in h:
        #                     d[k] = h[k][()]

        mask = mask if mask is not None else ()

        # TODO: Read only relevant fields
        d = {}
        with h5py.File(filename, 'r') as h:
            for k in h.keys():
                dd = h[k][mask]
                if k in self.column_mapping:
                    d[self.column_mapping[k]] = dd
                else:
                    d[k] = dd

        sim = self._create_catalog(name=name)
        sim._set_data(d)

        return sim
