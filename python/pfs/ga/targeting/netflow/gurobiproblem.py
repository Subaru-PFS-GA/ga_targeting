import os
import numpy as np

from .setup_logger import logger

try:
    import gurobipy as gbp
    from gurobipy import LinExpr
except ModuleNotFoundError:
    logger.warning('Module `gurobipy` is not available.')
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
            if isinstance(self.solver_options, dict):
                solver_options = self.solver_options
            else:
                solver_options = self.solver_options.to_dict()
                
            for key, value in solver_options.items():
                model.setParam(key, value)

        self.__cost = gbp.LinExpr()
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
        if lo is None:
            lo = -gbp.GRB.INFINITY

        if hi is None:
            hi = gbp.GRB.INFINITY

        if lo == 0 and hi == 1:
            vars = self.__model.addVars(indexes, name=name, vtype=gbp.GRB.BINARY)
        else:
            vars = self.__model.addVars(indexes, lb=lo, ub=hi, name=name, vtype=gbp.GRB.INTEGER)

        if cost is not None:
            self.add_cost(gbp.LinExpr(cost, [ vars[i] for i in indexes ]))
            
        self._variables[name] = vars
        return vars

    def get_variables(self):
        return self.__model.getVars()
    
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

    def get_constraints(self):
        return self.__model.getConstrs()

    def solve(self):
        self.__model.setObjective(self.__cost)
        self.__model.update()
        self.__model.optimize()

        if self.__model.Status == gbp.GRB.INFEASIBLE:
            logger.error('Model is infeasible which means certain constraints are not satisfiable. '
                         'Below is the list of constraints that are not satisfiable. Turn on the '
                         '`use_named_variables` config option to see the names of the constraints.')

            self.__model.computeIIS()
            constrs = self.__model.getConstrs()
            iisidx = np.where(self.__model.IISConstr)[0]
            for i in iisidx:
                self._infeasible_constraints[constrs[i].ConstrName] = constrs[i]
                logger.error(f'Infeasible constraint: `{constrs[i].ConstrName}`.')

        return self.__model.Status == gbp.GRB.OPTIMAL
    
    def update(self):
        self.__model.update()

    def write_problem(self, filename):
        # Assume that file extension is of a model file (such as mps),
        # otherwise, Gurobi will write something else to the file.
        self.__model.write(filename)

    def read_problem(self, filename):
        if os.path.isfile(filename):
            self.__model = gbp.read(filename)
        else:
            raise FileNotFoundError(f'File `{filename}` does not exist.')
        
        # TODO: populate variable list

    def write_solution(self, filename):
        self.__model.write(filename)

    def read_solution(self, filename):
        self.__model.read(filename)