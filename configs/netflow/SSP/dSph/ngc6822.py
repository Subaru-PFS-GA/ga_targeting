PROPOSALID = 'S25B-OT02'
CATID_SKY_GAIA = 1006
CATID_SKY_PS1 = 1007
CATID_FLUXSTD = 3006
CATID_SCIENCE_CO = 10091
CATID_SCIENCE_GA = 10092
CATID_SCIENCE_GE = 10093

extra_columns = {
    'proposalid': dict(
        constant = PROPOSALID,
        dtype = 'string',
    ),
    'obcode': dict(
        pattern = "SSP_GA_S25A_{field}_{target_list}_{{targetid:d}}_{resolution}",
        dtype = 'string'
    )
}

config = dict(
    targets = {
        "hsc": dict(
            path = "/datascope/subaru/user/dobos/netflow/catalogs/scl_obs.feather",
            # reader = None
            reader_args = dict(),
            column_map = {'objid': 'targetid'},
            prefix = "sci",
            frame='icrs',
            epoch = 2016.0,
            catid = CATID_SCIENCE_GA,
            extra_columns = extra_columns,
            photometry = dict(
                filters = {
                    "g_hsc": dict(
                        mag = 'obs_hsc_g',
                        mag_err = 'err_hsc_g',
                    ),
                    "i2_hsc": dict(
                        mag = 'obs_hsc_i2',
                        mag_err = 'err_hsc_i2',
                    ),
                    "nb515_hsc": dict(
                        mag = 'obs_hsc_nb515',
                        mag_err = 'err_hsc_nb515',
                    ),
                }
            ),
        ),
        # Miho
        "sky": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sextans/sky_sculptor.feather",
            reader_args = dict(),
            column_map = {
                'sky_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
            },
            prefix = "sky",
            catid = CATID_SKY_PS1,
            extra_columns = extra_columns,
        ),
        # Miho
        "fluxstd": dict(
            path = "$PFS_TARGETING_DATA/data/targeting/dSph/sextans/fluxstd_sculptor.feather",
            reader_args = dict(),
            column_map = {
                'fluxstd_id': 'targetid',
                'ra': 'RA',
                'dec': 'Dec',
                'parallax_error': 'err_parallax',
                'pmra_error': 'err_pmra',
                'pmdec_error': 'err_pmdec',
            },
            prefix = "cal",
            frame = 'icrs',
            epoch = 2016.0,
            catid = CATID_FLUXSTD,
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
                    'ps1_g': [16, 19.5],
                    'ps1_g-ps1_r': [0.25, 0.4],
                }
            )
        ),
    }
)
