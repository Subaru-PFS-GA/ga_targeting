import numpy as np

config = dict(
    # Munoz
    # obs_path = "$PFS_TARGETING_DATA/data/targeting/dSph/bootes/bootes_cfht.csv",
    # pmap_path = "$PFS_TARGETING_DATA/data/targeting/dSph/bootes/pmap/bootes",
    
    # Sato
    obs_path = "$PFS_TARGETING_DATA/data/targeting/dSph/bootes/bootes_hsc_sato_v2.csv",
    pmap_path = "$PFS_TARGETING_DATA/data/targeting/dSph/bootes/pmap/bootes_hsc",
    
    isochrones_path = "$PFS_TARGETING_DATA/data/cmdfit/isochrones/dartmouth/import/afep0_cfht_sdss_hsc",
    cut_nb = False,
    keep_blue = False,
    lp_member_limit = None,
    gaia_crossmatch = True,
    gaia_crossmatch_radius = 0.1
)