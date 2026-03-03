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
    extents = [[-0.1, 1.0], [14.0, 22.0]],
    bins = [100, 100],
    population_names = [
        'thin1', 'thin2', 'thin3', 'thick', 'halo', 'GC'
    ],
    population_weights = list(population_weights),
    merge_list = [list(range(10)), list(range(10, 12))],
)