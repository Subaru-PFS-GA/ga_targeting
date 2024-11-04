import os
from unittest import TestCase
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow.compat.v2 as tf

import  pfs.ga.targeting
from pfs.ga.targeting.selection.andselection import AndSelection
from pfs.ga.targeting.io import ObservationSerializer, Hdf5SimulationReader
from pfs.ga.targeting.data import Observation, Simulation
from pfs.ga.targeting.photometry import *
from pfs.ga.targeting.instrument import SubaruHSC
from pfs.ga.targeting.diagram import *
from pfs.ga.targeting.diagram.diagram import Diagram
from pfs.ga.targeting.diagram.xyaxis import XYAxis
from pfs.ga.targeting.projection import Pointing, WcsProjection
from pfs.ga.targeting.selection import *
from pfs.ga.targeting import Isochrone
from pfs.ga.isochrones.dartmouth import Dartmouth

class TestBase(TestCase):
    @classmethod
    def setUpClass(cls):
        tf.enable_v2_behavior()
        gpus = tf.config.list_physical_devices('GPU') 
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        try:
            tf.compat.v1.enable_eager_execution()
        except ValueError:
            pass

    def setUp(self):
        plt.figure(figsize=(10, 6))

        def get_env(k):
            return os.environ[k].strip('"') if k in os.environ else ''

        self.PFS_TARGETING_ROOT = get_env('PFS_TARGETING_ROOT')
        self.PFS_TARGETING_DATA = get_env('PFS_TARGETING_DATA')
        self.PFS_TARGETING_TEMP = get_env('PFS_TARGETING_TEMP')

    def get_test_data_file(self, filename):
        return os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data', filename)

    def get_filename(self, ext):
        filename = type(self).__name__[:-4] + '_' + self._testMethodName[5:] + ext
        return filename

    def save_fig(self, f=None, filename=None):
        if f is None:
            f = plt
        if filename is None:
            filename = self.get_filename('.png')

        f.tight_layout()
        f.savefig(os.path.join(self.PFS_TARGETING_TEMP, filename))

    def load_test_observation(self) -> Observation:
        fn =  os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi.feather')
        s = ObservationSerializer()
        return s.read(fn)

    def load_test_simulation(self) -> Simulation:
        fn = '/datascope/subaru/data/cmdfit/run/umi/sim/mix_bin_200k_hsc_007/sample.h5'
        r = Hdf5SimulationReader()
        cm = {}
        for prefix in ['', 'obs_', 'err_', 'flux_', 'obs_flux_', 'err_flux_', 'counts_', 'obs_counts_', 'err_counts_']:
            cm[prefix + 'hsc_g2'] = prefix + 'hsc_g'
            cm[prefix + 'hsc_i2'] = prefix + 'hsc_i'
        r.column_mapping = cm
        r.append_photometry(SubaruHSC.photometry())
        return r.read(fn)

    def load_test_isogrid(self):
        isogrid = Dartmouth()
        isogrid.load(os.path.join(self.PFS_TARGETING_DATA, 'isochrones/dartmouth/import/afep0_cfht_sdss_hsc_nb_bosz/isochrones.h5'))
        return isogrid

    def get_test_plot(self, projection=None):
        f, ax = plt.subplots(1, 1, figsize=(3.5, 3.5), dpi=250, subplot_kw=dict(projection=projection))
        return f, ax

    def get_test_diagram(self):
        d = Diagram([XYAxis('x'), XYAxis('y')])
        return d

    def get_test_cmd(self):
        photometry = SubaruHSC.photometry()
      
        cmd = CMD([
            ColorAxis(
                Color([photometry.magnitudes['g'], photometry.magnitudes['i']]),
                limits=(-1.5, 4.5)),
            MagnitudeAxis(
                photometry.magnitudes['g'],
                limits=(15.5, 24.5))
        ])

        return cmd, photometry

    def get_test_ccd(self):
        photometry = SubaruHSC.photometry()
      
        ccd = CCD([
            ColorAxis(
                Color([photometry.magnitudes['g'], photometry.magnitudes['i']]),
                limits=(-1.5, 4.5)),
            ColorAxis(
                Color([photometry.magnitudes['nb515'], photometry.magnitudes['g']]),
                limits=(-2, 2))
        ])

        return ccd, photometry

    def get_test_magnitude_selection(self, cmd):
        return MagnitudeSelection(cmd.axes[1], 18.5, 22.5)

    def get_test_color_selection(self, cmd):
        return ColorSelection(cmd.axes[0], 0, 1.5)
        
    def get_test_selection(self, cmd):
        return AndSelection([
            self.get_test_magnitude_selection(cmd),
            self.get_test_color_selection(cmd)
        ])

    def get_test_isochrone(self):
        isogrid = self.load_test_isogrid()
        photometry = SubaruHSC.photometry()

        name_mappings = {
            # 'hsc_g': 'hsc_g2',
            # 'hsc_i': 'hsc_i2'
        }

        iso = Isochrone()
        iso.from_isogrid(photometry, isogrid, Fe_H=-1.75, log_t=10.135, DM=19.2, name_mappings=name_mappings)

        return iso

    def get_projection(self, obs):
        ra, dec = obs.get_coords()
        p = WcsProjection(Pointing([ ra.mean(), dec.mean() ], obs_time=datetime(2020, 1, 1, 12, 0, 0)), proj='TAN')
        return p