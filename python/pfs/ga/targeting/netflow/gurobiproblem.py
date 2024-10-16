import logging
try:
    import gurobipy as gbp
    from gurobipy import LinExpr
except ModuleNotFoundError:
    logging.warning('Module `gurobipy` is not available.')
    gbp = None

from .ilpproblem import ILPProblem

class GurobiProblem(ILPProblem):
    def __init__(self, name='problem', solver_options=None):
        super().__init__(name=name, solver_options=solver_options)

        self.__model = None
        self.__cost = None
        self.__sum = None

        self.__create_model()

    def __get_cost(self):
        return self.__cost
    
    cost = property(__get_cost)

    def __get_sum(self):
        return self.__sum
    
    sum = property(__get_sum)

    def __create_model(self):
        model = gbp.Model(self.name)
        model.ModelSense = 1  # Minimize
        if self.solver_options is not None:
            for key, value in self.solver_options.items():
                model.setParam(key, value)

        self.__cost = model.addVar(name="cost", vtype=gbp.GRB.CONTINUOUS)
        self.__sum = gbp.quicksum
        self.__model = model

    def add_variable(self, name, lo, hi):       
        if lo is None:
            lo = -gbp.GRB.INFINITY
        if hi is None:
            hi = gbp.GRB.INFINITY
        if lo == 0 and hi == 1:
            var = self.__model.addVar(name=name, vtype=gbp.GRB.BINARY)
        else:
            var = self.__model.addVar(lb=lo, ub=hi, name=name, vtype=gbp.GRB.INTEGER)
        self._variables[name] = var
        return var
    
    def add_variable_array(self, name, indexes, lo, hi, cost=None):
        if cost is None:
            cost = 0.0

        if lo is None:
            lo = -gbp.GRB.INFINITY

        if hi is None:
            hi = gbp.GRB.INFINITY

        if lo == 0 and hi == 1:
            vars = self.__model.addVars(indexes, name=name, obj=cost, vtype=gbp.GRB.BINARY)
        else:
            vars = self.__model.addVars(indexes, lb=lo, ub=hi, name=name, obj=cost, vtype=gbp.GRB.INTEGER)
            
        self._variables[name] = vars
        return vars
    
    def get_value(self, variable):
        return variable.X
    
    def add_cost(self, additional_cost):
        self.__cost += additional_cost
    
    def add_constraint(self, name, constraint, lazy=None):
        if lazy is not None:
            constraint.Lazy = lazy
        self._constraints[name] = constraint

        # TODO: replace with addLConstr to speed up model building
        #       it requires rewriting everything because right now
        #       constraints are created as python expressions
        # self.__model.addLConstr(lhs, sense, rhs, name=name)

        self.__model.addConstr(constraint, name=name)

    def add_lazy_constraint(self, name, constraint):
        self.add_constraint(name, constraint, lazy=1)

    def add_linear_constraint(self, name, coeffs, variables, sense, rhs):
        self.__model.addLConstr(LinExpr(coeffs, variables), sense, rhs, name=name)

    def solve(self):
        self.__model.setObjective(self.__cost)
        self.__model.update()
        self.__model.optimize()
    
    def update(self):
        self.__model.update()

    def write_problem(self, filename):
        self.__model.write(filename)

    def read_problem(self, filename):
        self.__model.read(filename)
        # TODO: populate variable list

    def write_solution(self, filename):
        self.__model.write(filename)

    def read_solution(self, filename):
        self.__model.read(filename)