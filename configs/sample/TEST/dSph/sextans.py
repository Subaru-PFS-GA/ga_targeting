import numpy as np

config = dict(
    obs_path = "$PFS_TARGETING_DATA/data/targeting/dSph/sextans/sextans_tpall3e_g24.cat",
    pmap_path = "$PFS_TARGETING_DATA/data/targeting/dSph/sextans/pmap/TEST/sextans_TEST_003/",
    isochrones_path = "$PFS_TARGETING_DATA/data/cmdfit/isochrones/dartmouth/import/afep0_cfht_sdss_hsc",
    isochrones_name_mapping = {
        'hsc_r2': 'hsc_r',
        'hsc_i2': 'hsc_i',
    },
    cut_nb = True,
    keep_blue = True,
    lp_member_limit = np.log(0.001),
    gaia_crossmatch = True,
    gaia_crossmatch_radius = 0.1
)