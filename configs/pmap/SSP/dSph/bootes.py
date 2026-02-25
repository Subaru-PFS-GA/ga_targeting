import numpy as np

population_weights = np.array([
    0.0141608 , 0.00962137,
    0.02478971, 0.0164748 ,
    0.04102093, 0.02646154,
    0.25325579, 0.16523518,
    0.05996296, 0.03948173,
    0.209696  , 0.13983918])
population_weights /= population_weights.sum()

config = dict(
    cut_nb = False,
    keep_blue = False,
    extents = [[0.0, 1.5], [16.8, 23.0]],
    bins = [100, 100],
    population_names = [
        'thin1', 'thin2', 'thin3', 'thick', 'halo', 'dSph'    
    ],
    population_weights = list(population_weights),
    merge_list = [list(range(10)), list(range(10, 12))],
    sim_path = "/datascope/subaru/user/khayasi/cmdfit/run/bootes/sim/bin_chab_sdss_250k_004/",
    
    # sim_path = "/datascope/subaru/data/cmdfit/run/umi/sim/bin_chab_nb_uni_250k_001",
    # use_p_stars = True,
    
    # obs_path = "/datascope/subaru/data/catalogs/Munoz+18/munoz.h5"
    obs_path = "/datascope/subaru/data/targeting/dSph/bootes/bootes_cfht.csv"
    
    #obs_path = "/home/khayasi/Subaru-PFS-GA/ga_hscdata/sextans_within_2deg.csv"
)