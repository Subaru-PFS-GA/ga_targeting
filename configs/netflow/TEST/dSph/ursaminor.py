from datetime import datetime, timedelta
from pfs.ga.targeting.targets.dsph import Fornax
from pfs.ga.targeting.instrument import SubaruHSC

extra_columns = {
    'proposalid': dict(
        pattern = "SSP_GA_{obs_time:%Y%m%d}_{name}",
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    pointings = [
        dict(ra=229.2, dec=67.90, posang=0),
    ],
    targets = {
        "hsc": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/ursaminor_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            frame='icrs',
            catid = 10088,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_hsc": dict(
                        mag = 'obs_hsc_g',
                        mag_err = 'err_hsc_g',
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
            )
        ),
        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/sky_ursaminor.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = 1007,
            extra_columns = extra_columns,
        ),
        "fluxstd": dict(
            # path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/fluxstd_ursaminor.feather",
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/PS1_UMI_fluxstd_2_dobos.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                # 'ra': 'RA',
                # 'dec': 'Dec',
                'parallax_error': 'err_parallax',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
                'radial_velocity': 'rv',
                'radial_velocity_error': 'err_rv',
            },
            # mask = 'lambda df: df["prob_f_star"] > 0.5',
            # mask = 'lambda df: df["prob_f_star"] > 0.1',
            prefix = "cal",
            frame = 'icrs',
            epoch = 2016.0,
            catid = 3006,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    'g_ps1': dict(
                        mag = 'gPSFMag',
                        mag_err = 'gPSFMagErr'
                    ),
                    'r_ps1': dict(
                        mag = 'rPSFMag',
                        mag_err = 'rPSFMagErr'
                    ),
                    'i_ps1': dict(
                        mag = 'iPSFMag',
                        mag_err = 'iPSFMagErr'
                    ),
                    'z_ps1': dict(
                        mag = 'zPSFMag',
                        mag_err = 'zPSFMagErr'
                    ),
                    'y_ps1': dict(
                        mag = 'yPSFMag',
                        mag_err = 'yPSFMagErr'
                    )
                }
                # bands = {
                #     b: dict(
                #         filter = f'filter_{b}',                     # Column storing filter names
                #         # psf_mag = f'psf_mag_{b}',                 # All None in data file
                #         # psf_mag_err = f'psf_mag_error_{b}',       # All None in data file
                #         psf_flux = f'psf_flux_{b}',
                #         psf_flux_err = f'psf_flux_error_{b}',
                #     ) for b in 'gri'
                # },
                # limits = {
                #     'ps1_g': [13, 18],
                # }
            )
        ),
        "guide": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/guide_ursaminor.feather",
            reader_args = dict(),
            prefix = "ag",
            photometry = dict(
                filters = {
                    "gaia_g": dict(
                        flux = "flux_gaia_g",
                        flux_err = "err_flux_gaia_g",
                    ),
                    "gaia_bp": dict(
                        flux = "flux_gaia_bp",
                        flux_err = "err_flux_gaia_bp",
                    ),
                    "gaia_rp": dict(
                        flux = "flux_gaia_rp",
                        flux_err = "err_flux_gaia_rp",
                    )
                }
            )
        )
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
