from datetime import datetime, timedelta
from pfs.ga.targeting.targets.dsph import Fornax
from pfs.ga.targeting.instrument import SubaruHSC

config = dict(
    targets = {
        "dsph": dict(
            path = "/datascope/subaru/data/targeting/dSph/ursaminor/ursaminor_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            epoch = "J2016.0",
            catid = 15001,
            proposalid = "PFS-2024-001",
            filters = {
                "g_hsc": dict(
                    mag = 'obs_hsc_g',
                    mag_err = 'err_hsc_g',
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
                    mag = 'obs_hsc_i',
                    mag_err = 'err_hsc_i',
                ),
                "nb515_hsc": dict(
                    mag = 'obs_hsc_nb515',
                    mag_err = 'err_hsc_nb515',
                ),
            }
        ),
        "sky": dict(
            path = "/datascope/subaru/data/targeting/dSph/ursaminor/sky_ursaminor.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky"
        ),
        "fluxstd": dict(
            path = "/datascope/subaru/data/targeting/dSph/ursaminor/fluxstd_ursaminor.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "cal",
            bands = {
                b: dict(
                    filter = f'filter_{b}',                   # Column storing filter names
                    psf_mag = f'psf_mag_{b}',
                    psf_mag_err = f'psf_mag_error_{b}',
                    psf_flux = f'psf_flux_{b}',
                    psf_flux_err = f'psf_flux_error_{b}',
                ) for b in 'gri'
            }
        ),
    },
    # Override the minimum number of calibration targets
    netflow_options = dict(
        cobra_groups = {
            'cal_location': dict(
                min_targets = 0,
            )
        }
    )
)
