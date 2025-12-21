import numpy as np

config = dict(
    cut_nb = True,
    keep_blue = True,
    extents = [[0.1, 2.0], [17.0, 23.5]],
    bins = [100, 100],
    population_names = [
        'thin1', 'thin2', 'thin3', 'thick', 'halo', 'dSph'    
    ],
    population_weights = [
        2.02985435e-04, 8.36941284e-05,
        3.63036869e-02, 1.50133136e-02,
        5.82054620e-04, 2.46598771e-04,
        5.40316250e-02, 2.28643935e-02,
        1.33522569e-01, 5.75283004e-02,
        4.74330132e-01, 2.05290646e-01
    ],
    merge_list = [list(range(10)), list(range(10, 12))],
    sim_path = "/datascope/subaru/user/khayasi/cmdfit/run/sext/sim/bin_chab_nb_250k_001/",
    
    # sim_path = "/datascope/subaru/data/cmdfit/run/umi/sim/bin_chab_nb_uni_250k_001",
    # use_p_stars = True,
    
    obs_path = "$PFS_TARGETING_DATA/data/targeting/dSph/sextans/sextans_tpall3e_g24.cat"
    # obs_path = "/home/khayasi/Subaru-PFS-GA/ga_hscdata/sextans_within_1deg.csv"
)