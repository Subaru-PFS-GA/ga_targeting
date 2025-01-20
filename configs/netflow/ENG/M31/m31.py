from datetime import datetime, timedelta

config = dict(
    targets = {
        "m31": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/M31/M31_ENG/m31_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            # epoch = "J2000.0",
            catid = 10088,
            extra_columns = {
                'proposalid': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}",
                    dtype = 'string',
                ),
                'obcode': dict(
                    pattern = "SSP_GA_ENG_{obs_time:%Y%m%d}_{name}_{{targetid:d}}_{resolution}",
                    dtype = 'string'
                )
            },
            photometry = dict(
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
        ),
        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/M31/M31_ENG/sky_m31.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = 2007,
        ),
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/M31/M31_ENG/fluxstd_m31.feather",
            reader_args = dict(),
            column_map = {
                'obj_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "cal",
            catid = 3006,
            photometry = dict(
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                   # Column storing filter names
                        # psf_mag = f'psf_mag_{b}',
                        # psf_mag_err = f'psf_mag_error_{b}',
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'gri'
                }
            ),
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
