from ..config import Config

class GurobiOptionsConfig(Config):
    def __init__(self):
        # Random seed
        self.seed = None

        # Agressiveness of presolve which tries to eliminate variables from the LP problem
        # -1 means automatic, which is the best choice to prevent the presolver from
        # taking too much time before the heuristics kick in
        self.presolve = -1
        
        # 3 means concurrent, 4 means deterministic concurrent. Set it to 4 is 
        # reproducibility is needed
        self.method = 3

        # Degenerate simplex moves, set to 0 to prevent too much time to be spent on
        # trying to improve the current solution
        self.degenmoves = 0

        # How much of the time to spend by performing heuristics
        self.heuristics = 0.5
        
        # mipfocus=1 is balanced toward finding more feasible solutions
        # mipfocus=2 is balanced toward proving that the current solution is the best
        # mipfocus=3 is to be used when the objection bound is moving very slowly
        self.mipfocus = 1

        # Relative stopping criterion for bounds on the objective. Since it is expressed
        # relative to the objective value, it should be chosen accordingly. Run the
        # optimization to see what the typical objective value is and set the value of
        # mipgap to be in the range of the lowest priority science targets.
        self.mipgap = 0.01

        self.mipgapabs = None

        # The number of threads to use for the optimization. Set to 0 to determine the
        # thread count automatically.
        self.threads = 0

        self.LogToConsole = 1
        
        self.timelimit = 300

        # When set to 1, force checking for integrality of integer and binary variables
        self.IntegralityFocus = 0

        # Tolerance for integrality of integer and binary variables
        self.IntFeasTol = 1e-5

        super().__init__()