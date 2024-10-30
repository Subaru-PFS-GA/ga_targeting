class ILPProblem():
    def __init__(self, name='problem', solver_options=None):
        self.__name = name
        self.__solver_options = solver_options if solver_options is not None else {}
        self._variables = {}
        self._constraints = {}
        self._infeasible_constraints = {}

    def __get_name(self):
        return self.__name
    
    name = property(__get_name)

    def __get_solver_options(self):
        return self.__solver_options
    
    solver_options = property(__get_solver_options)

    def __get_variables(self):
        return self._variables
    
    variables = property(__get_variables)
    
    def __get_constraints(self):
        return self._constraints
    
    constraints = property(__get_constraints)

    def __get_infeasible_constraints(self):
        return self._infeasible_constraints
    
    infeasible_constraints = property(__get_infeasible_constraints)

    def add_variable(self, name, lo, hi):
        raise NotImplementedError()
    
    def add_constraint(self, name, constraint):
        raise NotImplementedError()
    
    def add_lazy_constraint(self, name, constraint):
        raise NotImplementedError()
    
    def add_linear_constraint(self, name, coeffs, variables, sense, rhs):
        raise NotImplementedError()
    
    def solve(self):
        raise NotImplementedError()
    
    def update(self):
        raise NotImplementedError()
    
    def write_problem(self, filename):
        raise NotImplementedError()

    def read_problem(self, filename):
        raise NotImplementedError()

    def write_solution(self, filename):
        raise NotImplementedError()

    def read_solution(self, filename):
        raise NotImplementedError()
