from ..config import Config

class PMapConfig(Config):
    def __init__(self,
                 cut_nb: bool = False,
                 keep_blue: bool = False,
                 extents = [[0.1, 2.0], [17.0, 23.5]],
                 bins = [100, 100],
                 population_names = None,
                 population_weights = None,
                 merge_list = None):
        
        # Use the NB515 filter to define the color cuts
        self.cut_nb = cut_nb

        # TODO
        self.keep_blue = keep_blue

        # Probability map extents
        self.extents = extents

        # Histogram bins for the CMD
        self.bins = bins

        # Whether binary stars are simulated separately
        self.binaries = True

        # Population names
        self.population_names = []

        # Override weights for each population
        self.population_weights = population_weights

        # Population merge list
        self.merge_list = merge_list

        # Simulation path where sample.h5 and sim.h5 are located
        self.sim_path = None

        # Observations file, only for generating the reports
        self.obs_path = None

        super().__init__()