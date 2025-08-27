import numpy as np

config = dict(
    cut_nb = True,
    keep_blue = False,
    extents = [[0.1, 3.0], [17.0, 24.5]],
    bins = [100, 100],
    population_names = [
        'thin1', 'thin2', 'thin3', 'thick', 'halo', 'm31'
    ],
    population_weights = [
        2.02985435e-09, #8.36941284e-05,
        3.63036869e-09, #1.50133136e-02,
        5.82054620e-09, #2.46598771e-04,
        5.40316250e-09, #2.28643935e-02,
        2.33522569e-01, #5.75283004e-02,
        1.74330132e-00, #2.05290646e-00
    ],
    merge_list = [list(range(5)), list(range(5, 6))],
    sim_path = "/datascope/subaru/data/cmdfit/run/m31/sim/nobin_chab_250k_002",
    
    # sim_path = "/datascope/subaru/data/cmdfit/run/umi/sim/bin_chab_nb_uni_250k_001",
    # use_p_stars = True,
    
    obs_path = "/datascope/subaru/data/targeting/m31/M31Catalog_forPFS.csv"
)

config['population_weights'] /= np.sum(config['population_weights'])