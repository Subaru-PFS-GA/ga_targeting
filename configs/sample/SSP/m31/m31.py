import numpy as np

config = dict(
    obs_path = "$PFS_TARGETING_DATA/data/targeting/m31/m31_hsc_catalog_for_pfs_mlprob.csv",
    pmap_path = "$PFS_TARGETING_DATA/data/targeting/m31/m31_sector_0_SSP/pmap/m31_sector_0_SSP_001",
    isochrones_path = "$CMDFIT_DATA/isochrones/dartmouth/import/afep0_cfht_sdss_hsc",
    cut_nb = True,
    keep_blue = False,
    lp_member_limit = 0,
    gaia_crossmatch = True,
    gaia_crossmatch_radius = 0.1
)