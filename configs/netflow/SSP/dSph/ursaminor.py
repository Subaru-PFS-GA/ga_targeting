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
    # pointings = [
    #     dict(ra=229.2, dec=67.90, posang=0),
    # ],
    # Override the minimum number of calibration targets
    netflow_options = dict(
        target_classes = {
            'sky': dict(
                prefix = 'sky',
                min_targets = 240,
                max_targets = 320,
            ),
            'cal': dict(
                prefix = 'cal',
                min_targets = 40,
                max_targets = 240,
            ),
        },
        cobra_groups = {
            'cal_location': dict(
                min_targets = 0,
            )
        }
    ),
    targets = {
        "hsc": dict(
            # path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/ursaminor_obs.feather",
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/priority/ursaminor_nb_anc/hsc_umi_priorities.feather",
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
        "anc": dict(
            # path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/ursaminor_obs.feather",
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/umi.anc.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            frame='icrs',
            catid = 10088,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_ps": dict(
                        mag = 'obs_ps_g',
                        mag_err = 'err_ps_g',
                    ),
                    "r_ps": dict(
                        mag = 'obs_ps_r',
                        mag_err = 'err_ps_r',
                    ),
                    "i_ps": dict(
                        mag = 'obs_ps_i',
                        mag_err = 'err_ps_i',
                    ),
                    "z_ps": dict(
                        mag = 'obs_ps_z',
                        mag_err = 'err_ps_z',
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

        # DOBOS
        # "fluxstd": dict(
        #     path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/PS1_UMI_fluxstd_2_dobos.feather",
        #     reader_args = dict(),
        #     column_map = {
        #         'obj_id': 'targetid',
        #         # 'ra': 'RA',
        #         # 'dec': 'Dec',
        #         'parallax_error': 'err_parallax',
        #         'pmra_error': 'err_pmra',
        #         'pmdec_error': 'err_pmdec',
        #         'radial_velocity': 'rv',
        #         'radial_velocity_error': 'err_rv',
        #     },
        #     # mask = 'lambda df: df["prob_f_star"] > 0.5',
        #     # mask = 'lambda df: df["prob_f_star"] > 0.1',
        #     prefix = "cal",
        #     frame = 'icrs',
        #     epoch = 2016.0,
        #     catid = 3006,
        #     extra_columns = extra_columns,
        #     photometry = dict(
        #         filters = {
        #             'g_ps1': dict(
        #                 mag = 'gPSFMag',
        #                 mag_err = 'gPSFMagErr'
        #             ),
        #             'r_ps1': dict(
        #                 mag = 'rPSFMag',
        #                 mag_err = 'rPSFMagErr'
        #             ),
        #             'i_ps1': dict(
        #                 mag = 'iPSFMag',
        #                 mag_err = 'iPSFMagErr'
        #             ),
        #             'z_ps1': dict(
        #                 mag = 'zPSFMag',
        #                 mag_err = 'zPSFMagErr'
        #             ),
        #             'y_ps1': dict(
        #                 mag = 'yPSFMag',
        #                 mag_err = 'yPSFMagErr'
        #             )
        #         }
        #     )
        # ),

        # MIHO NEW
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/fluxstd_ursaminor_miho_20250228.feather",
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax_error': 'err_parallax',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
            },
            mask = 'lambda df: df["prob_f_star"] > 0.5',
            prefix = "cal",
            frame = 'icrs',
            epoch = 2016.0,
            catid = 3006,
            extra_columns = extra_columns,
            photometry = dict(
                bands = {
                    b: dict(
                        filter = f'filter_{b}',                     # Column storing filter names
                        # psf_mag = f'psf_mag_{b}',                 # All None in data file
                        # psf_mag_err = f'psf_mag_error_{b}',       # All None in data file
                        psf_flux = f'psf_flux_{b}',
                        psf_flux_err = f'psf_flux_error_{b}',
                    ) for b in 'grizy'
                },
                limits = {
                    'ps1_g': [16, 19],
                    'ps1_g-ps1_r': [0.25, 0.4],
                }
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
)
