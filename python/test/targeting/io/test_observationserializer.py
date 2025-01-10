import os

from test_base import TestBase
import pfs.ga.targeting
from pfs.ga.targeting.config import TargetListConfig
from pfs.ga.targeting.instrument.subaruhsc import SubaruHSC
from pfs.ga.targeting.io import ObservationSerializer

class ObservationSerializerTest(TestBase):
    def test_read_hdf5(self):
        fn = '/datascope/subaru/data/cmdfit/catalog/Munoz+18/munoz.h5'
        r = ObservationSerializer()
        r.read(fn, 'obs/umi/cfht')

    def test_read_hsc(self):
        fn = '/datascope/subaru/data/cmdfit/dSph/ursaminor_tpall3e_g24.cat'
        r = SubaruHSC.text_observation_reader()
        obs = r.read(fn)

    def test_read_with_filters(self):
        config = dict(
            path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi.feather'),
            reader_args = dict(),
            column_map = {
                'objid': 'targetid'
            },
            prefix = "sci",
            epoch = "J2000.0",
            catid = 15001,
            proposalid = "Test",
            filters = {
                "g_hsc": dict(
                    mag = 'sdss_g',
                    mag_err = 'err_sdss_g',
                    # flux = "g_hsc_flux",
                    # flux_err = "g_hsc_flux_err",
                    # psf_mag = "g_hsc",
                    # psf_mag_err = "g_hsc_err",
                    # psf_flux = "g_hsc_flux",
                    # psf_flux_err = "g_hsc_flux_err",
                    # fiber_mag = "g_hsc",
                    # fiber_mag_err = "g_hsc_err",
                    # fiber_flux = "g_hsc_flux",
                    # fiber_flux_err = "g_hsc_flux_err",
                    # total_mag = "g_hsc",
                    # total_mag_err = "g_hsc_err",
                    # total_flux = "g_hsc_flux",
                    # total_flux_err = "g_hsc_flux_err",
                ),
                "i_hsc": dict(
                    mag = 'sdss_r',
                    mag_err = 'err_sdss_r',
                ),
            }
        )

        target_list_config = TargetListConfig.from_dict(config)

        reader = ObservationSerializer(
            columns = target_list_config.columns,
            column_map = target_list_config.column_map,
            data_types = target_list_config.data_types,
            index = target_list_config.index,
            filters = target_list_config.filters,
            bands = target_list_config.bands,
            kwargs = target_list_config.reader_args,
        )
        catalog = reader.read(target_list_config.path)

        self.assertTrue('hsc' in catalog.photometry)
        self.assertTrue('g' in catalog.photometry['hsc'].magnitudes)
        self.assertTrue('i' in catalog.photometry['hsc'].magnitudes)

    def test_read_with_bands(self):
        config = dict(
            path = os.path.join(os.path.dirname(pfs.ga.targeting.__file__), '../../../../data/test/umi_fluxstd.feather'),
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'obj_id': 'orig_objid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax': 'parallax',
                'parallax_error': 'err_parallax',
                'pmra': 'pmra',
                'pmra_error': 'err_pmra',
                'pmdec': 'pmdec',
                'pmdec_error': 'err_pmdec',
            },
            prefix = "cal",
            bands = {
                b: dict(
                    filter = f'filter_{b}',                   # Column storing filter names
                    # psf_mag = f'psf_mag_{b}',
                    # psf_mag_err = f'psf_mag_error_{b}',
                    psf_flux = f'psf_flux_{b}',
                    psf_flux_err = f'psf_flux_error_{b}',
                    # fiber_mag = None,
                    # fiber_mag_err = None,
                    # fiber_flux = None,
                    # fiber_flux_err = None,
                    # total_mag = None,
                    # total_mag_err = None,
                    # total_flux = None,
                    # total_flux_err = None,
                ) for b in 'gr'
            }
        )
        
        target_list_config = TargetListConfig.from_dict(config)

        reader = ObservationSerializer(
            columns = target_list_config.columns,
            column_map = target_list_config.column_map,
            data_types = target_list_config.data_types,
            index = target_list_config.index,
            filters = target_list_config.filters,
            bands = target_list_config.bands,
            kwargs = target_list_config.reader_args,
        )
        catalog = reader.read(target_list_config.path)

        self.assertTrue('ps1' in catalog.photometry)
        self.assertTrue('g' in catalog.photometry['ps1'].magnitudes)
        self.assertTrue('r' in catalog.photometry['ps1'].magnitudes)