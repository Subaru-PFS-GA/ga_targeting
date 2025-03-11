import numpy as np

config = dict(
    obs_path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/ursaminor_tpall3e_g24.cat",
    pmap_path = "$PFS_TARGETING_DATA/data/targeting/dSph/ursaminor/pmap/ursaminor_nb",
    isochrones_path = "$PFS_TARGETING_DATA/data/cmdfit/isochrones/dartmouth/import/afep0_cfht_sdss_hsc",
    cut_nb = True,
    keep_blue = True,
    lp_member_limit = np.log(0.001),
    gaia_crossmatch = True,
    gaia_crossmatch_radius = 0.1
)