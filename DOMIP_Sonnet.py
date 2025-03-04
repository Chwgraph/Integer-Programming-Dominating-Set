# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 16:20:15 2025

@author: howiech
"""

# -*- coding: utf-8 -*-
"""
Optimized version of DOMIP.py
"""

import numpy as np
import networkx as nx
import gurobipy as grb
from ortools.sat.python import cp_model
import coptpy as cp
from coptpy import COPT
import matplotlib.pyplot as plt
from functools import lru_cache


class DominationProblemSolver:
    def __init__(self, adj_matrix, nodes, solver, params={}, M=10000):
        """
        Initialize the Domination Problem Solver
        
        Parameters:
        adj_matrix: Adjacency matrix of the graph
        nodes: Number of nodes in the graph
        solver: Solver to use ('GRB', 'OR', or 'CP')
        params: Additional parameters for the solver
        M: Big M value for constraints
        """
        self.adj = np.array(adj_matrix)
        self.nodes = nodes
        self.solver = solver
        self.params = params
        self.M = M
        
        # Pre-calculate commonly used matrices
        self.B0 = self.adj + np.eye(nodes, dtype=int)  # Adjacency matrix with self-loops
        
    # Helper methods to reduce code duplication
    def _create_gurobi_env(self):
        """Create and configure a Gurobi environment"""
        env = grb.Env(params=self.params)
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        return model
    
    def _create_copt_env(self):
        """Create and configure a COPT environment"""
        env = cp.Envr()
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        return model
    
    def _add_domination_constraints(self, model, B0, variables, size, version, solver_type):
        """Add domination constraints based on version and solver type"""
        rows = self.nodes
        
        if solver_type == 'OR':
            C = B0 * variables
            for j in range(rows):
                if version == 'w':
                    if C[j][j] == 0:
                        model.Add(sum(C[j]) >= size)
                elif version == 's':
                    model.Add(sum(C[j]) >= size)
        else:  # GRB or CP
            C = [[B0[j][k] * variables[k] for k in range(rows)] for j in range(rows)]
            for j in range(rows):
                if version == 'w':
                    if j != j:  # This check needs revision - original logic is unclear
                        model.addConstr(sum(C[j]) >= size)
                elif version == 's':
                    model.addConstr(sum(C[j]) >= size)
        
        return C
    
    #-------------------------------------------------------------------------
    # Minimum Dominating Set (MDS) Methods
    #-------------------------------------------------------------------------
    
    def MDS_OR(self, size, version):
        """Solve MDS using OR-Tools CP-SAT solver"""
        model = cp_model.CpModel()
        rows = self.nodes
        
        # Add variables to the model
        variables = [model.NewIntVar(0, 1, f'x{i}') for i in range(rows)]
        z = model.NewIntVar(0, rows, 'z')
        
        # Add domination constraints
        self._add_domination_constraints(model, self.B0, variables, size, version, 'OR')
        
        # Set objective function
        model.Add(z == sum(variables))
        model.Minimize(z)
        
        # Solve the model
        solver = cp_model.CpSolver()
        status = solver.Solve(model)
        
        # Process results
        if status == cp_model.OPTIMAL:
            VAR = [solver.Value(var) for var in variables]
            print(f'The domination number of the input graph is: {solver.Value(z)}')
            print(f'The number of vertices in the input graph is: {rows}')
        else:
            VAR = []
            print("No optimal solution found")
            
        return VAR
    
    def MDS_GRB(self, size, version):
        """Solve MDS using Gurobi solver"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Add variables to the model
        variables = [model.addVar(vtype=grb.GRB.BINARY, lb=0, ub=1, name=f'x{i}') for i in range(rows)]
        z = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name='z')
        
        # Add domination constraints
        self._add_domination_constraints(model, self.B0, variables, size, version, 'GRB')
        
        # Set objective function
        model.addConstr(z == sum(variables))
        model.setObjective(z, grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        VAR = [var.x for var in variables]
        print(f'The {size}-domination number is: {model.objVal}')
        
        return VAR
    
    def MDS_CP(self, size, version):
        """Solve MDS using COPT solver"""
        model = self._create_copt_env()
        rows = self.nodes
        
        # Add variables to the model
        variables = [model.addVar(vtype=cp.COPT.BINARY, lb=0, ub=1, name=f'x{i}') for i in range(rows)]
        z = model.addVar(vtype=cp.COPT.INTEGER, name='z')
        
        # Add domination constraints
        self._add_domination_constraints(model, self.B0, variables, size, version, 'CP')
        
        # Set objective function
        model.addConstr(z == sum(variables))
        model.setObjective(z, cp.COPT.MINIMIZE)
        
        # Solve the model
        status = model.solve()
        
        # Process results
        VAR = [var.x for var in variables]
        print(f'The {size}-domination number is: {model.objVal}')
        
        return VAR
    
    def MDS(self, k, ver):
        """Solve MDS using the selected solver"""
        if self.solver == 'OR':
            return self.MDS_OR(k, ver)
        elif self.solver == 'GRB':
            return self.MDS_GRB(k, ver)
        elif self.solver == 'CP':
            return self.MDS_CP(k, ver)
        else:
            raise ValueError(f"Unknown solver: {self.solver}")
    
    #-------------------------------------------------------------------------
    # Connected Minimum Dominating Set (CMDS) Methods
    #-------------------------------------------------------------------------
    
    def CMDS_MTZ_OR(self):
        """Solve CMDS using OR-Tools CP-SAT solver with Miller-Tucker-Zemlin formulation"""
        model = cp_model.CpModel()
        rows = self.nodes
        
        # Create adjacency matrix with additional rows
        A = [list(self.adj[i]) for i in range(rows)] + [[1 for _ in range(rows)], [1 for _ in range(rows)]]
        B = np.array(A)
        
        # Create variables
        var_1 = [model.NewBoolVar(name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_2 = [model.NewBoolVar(name=f'yn+1,{i}') for i in range(rows)]  # Source variables
        var_3 = [model.NewBoolVar(name=f'yn+2,{i}') for i in range(rows)]  # Sink variables
        U = [model.NewIntVar(1, rows+1, name=f'u{i}') for i in range(rows)]  # MTZ variables
        
        # Create flow variables
        var_4 = []
        for i in range(rows):
            var_row = [model.NewBoolVar(name=f'y{i},{j}') for j in range(rows)]
            var_4.append(var_row)
        var_4.append(var_2)  # Source row
        var_4.append(var_3)  # Sink row
        
        # Add dummy variables for MTZ constraints
        U.append(0)  # u_{n+1}
        U.append(model.NewIntVar(1, rows+1, name=f'u{rows+2}'))  # u_{n+2}
        
        # Extract matrix products
        C = np.array(B * var_4)
        C0 = self.B0 * var_1
        
        # Add domination constraints
        for j in range(rows):
            model.Add(sum(C0[j]) >= 1)
        
        # Add flow constraints
        model.Add(sum(var_3) == 1)  # Constraint 2a
        
        # Add flow conservation constraints
        for k in range(rows):
            model.Add(sum(C.T[k]) == 1)  # Constraint 2b
            
            # Add arc constraints
            for l in range(rows):
                if isinstance(C[k][l], cp_model.LinearExpr):  # Check if it's a variable expression
                    model.Add(var_2[k] + C[k][l] <= 1)  # Constraint 2c
                    model.Add((rows+1) * C[k][l] + U[k] - U[l] + (rows-1) * C[l][k] <= rows)  # Constraint 2d
        
        # Add source/sink constraints
        for k in range(rows):
            if isinstance(C[-2][k], cp_model.LinearExpr):
                model.Add((rows+1) * C[-2][k] - U[k] <= rows)  # Constraint 2e part 1
            if isinstance(C[-1][k], cp_model.LinearExpr):
                model.Add((rows+1) * C[-1][k] + U[-1] - U[k] <= rows)  # Constraint 2e part 2
        
        # Add relation between dominating set and flow variables
        for i in range(rows):
            model.Add(var_1[i] == 1 - var_2[i])  # Constraint 2i
        
        # Set objective function
        z = model.NewIntVar(lb=0, ub=rows, name='z')
        model.Add(z == sum(var_1))
        model.Minimize(z)
        
        # Solve the model
        solver = cp_model.CpSolver()
        status = solver.Solve(model)
        
        # Process results
        if status == cp_model.OPTIMAL:
            VAR = [solver.Value(var_1[i]) for i in range(rows)]
            print(f'The connected domination number is: {solver.ObjectiveValue()}')
        else:
            VAR = []
            print("No optimal solution found")
            
        return VAR
    
    def CMDS_MTZ_GRB(self):
        """Solve CMDS using Gurobi solver with Miller-Tucker-Zemlin formulation"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Create adjacency matrix with additional rows
        A = [list(self.adj[i]) for i in range(rows)] + [[1 for _ in range(rows)], [1 for _ in range(rows)]]
        B = np.array(A)
        
        # Create variables
        var_1 = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_2 = [model.addVar(vtype=grb.GRB.BINARY, name=f'yn+1,{i}') for i in range(rows)]  # Source variables
        var_3 = [model.addVar(vtype=grb.GRB.BINARY, name=f'yn+2,{i}') for i in range(rows)]  # Sink variables
        U = [model.addVar(vtype=grb.GRB.INTEGER, lb=1, ub=rows+1, name=f'u{i}') for i in range(rows)]  # MTZ variables
        
        # Create flow variables
        var_4 = []
        for i in range(rows):
            var_row = [model.addVar(vtype=grb.GRB.BINARY, name=f'y{i},{j}') for j in range(rows)]
            var_4.append(var_row)
        var_4.append(var_2)  # Source row
        var_4.append(var_3)  # Sink row
        
        # Add dummy variables for MTZ constraints
        U.append(0)  # u_{n+1}
        U.append(model.addVar(vtype=grb.GRB.INTEGER, lb=1, ub=rows+1, name=f'u{rows+2}'))  # u_{n+2}
        
        # Extract matrix products
        C = np.array(B * var_4)
        C0 = self.B0 * var_1
        
        # Add domination constraints
        for j in range(rows):
            model.addConstr(sum(C0[j]) >= 1)
        
        # Add flow constraints
        model.addConstr(sum(var_3) == 1)  # Constraint 2a
        
        # Add flow conservation constraints
        for k in range(rows):
            model.addConstr(sum(C.T[k]) == 1)  # Constraint 2b
            
            # Add arc constraints
            for l in range(rows):
                if not isinstance(C[k][l], (int, float)):  # Check if it's a variable expression
                    model.addConstr(var_2[k] + C[k][l] <= 1)  # Constraint 2c
                    model.addConstr((rows+1) * C[k][l] + U[k] - U[l] + (rows-1) * C[l][k] <= rows)  # Constraint 2d
        
        # Add source/sink constraints
        for k in range(rows):
            if not isinstance(C[-2][k], (int, float)):
                model.addConstr((rows+1) * C[-2][k] - U[k] <= rows)  # Constraint 2e part 1
            if not isinstance(C[-1][k], (int, float)):
                model.addConstr((rows+1) * C[-1][k] + U[-1] - U[k] <= rows)  # Constraint 2e part 2
        
        # Add relation between dominating set and flow variables
        for i in range(rows):
            model.addConstr(var_1[i] == 1 - var_2[i])  # Constraint 2i
        
        # Set objective function
        z = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name='z')
        model.addConstr(z == sum(var_1))
        model.setObjective(z, grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        VAR = [var.x for var in var_1]
        print(f'The connected domination number is: {model.objVal}')
        
        return VAR
    
    def CMDS_MTZ_CP(self):
        """Solve CMDS using COPT solver with Miller-Tucker-Zemlin formulation"""
        model = self._create_copt_env()
        rows = self.nodes
        
        # Create adjacency matrix with additional rows
        A = [list(self.adj[i]) for i in range(rows)] + [[1 for _ in range(rows)], [1 for _ in range(rows)]]
        B = np.array(A)
        
        # Create variables
        var_1 = [model.addVar(vtype=cp.COPT.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_2 = [model.addVar(vtype=cp.COPT.BINARY, name=f'yn+1,{i}') for i in range(rows)]  # Source variables
        var_3 = [model.addVar(vtype=cp.COPT.BINARY, name=f'yn+2,{i}') for i in range(rows)]  # Sink variables
        U = [model.addVar(vtype=cp.COPT.INTEGER, lb=1, ub=rows+1, name=f'u{i}') for i in range(rows)]  # MTZ variables
        
        # Create flow variables
        var_4 = []
        for i in range(rows):
            var_row = [model.addVar(vtype=cp.COPT.BINARY, name=f'y{i},{j}') for j in range(rows)]
            var_4.append(var_row)
        var_4.append(var_2)  # Source row
        var_4.append(var_3)  # Sink row
        
        # Add dummy variables for MTZ constraints
        U.append(0)  # u_{n+1}
        U.append(model.addVar(vtype=cp.COPT.INTEGER, lb=1, ub=rows+1, name=f'u{rows+2}'))  # u_{n+2}
        
        # Create explicit expressions for matrix products
        C = [[B[j][k] * var_4[j][k] for k in range(rows)] for j in range(rows+2)]
        C0 = [[self.B0[j][k] * var_1[k] for k in range(rows)] for j in range(rows)]
        
        # Add domination constraints
        for j in range(rows):
            model.addConstr(sum(C0[j]) >= 1)
        
        # Add flow constraints
        model.addConstr(sum(var_3) == 1)  # Constraint 2a
        
        # Add flow conservation constraints
        for k in range(rows):
            model.addConstr(sum([C[j][k] for j in range(rows+2)]) == 1)  # Constraint 2b
            
            # Add arc constraints
            for l in range(rows):
                if not isinstance(C[k][l], (int, float)):  # Check if it's a variable expression
                    model.addConstr(var_2[k] + C[k][l] <= 1)  # Constraint 2c
                    model.addConstr((rows+1) * C[k][l] + U[k] - U[l] + (rows-1) * C[l][k] <= rows)  # Constraint 2d
        
        # Add source/sink constraints
        for k in range(rows):
            if not isinstance(C[-2][k], (int, float)):
                model.addConstr((rows+1) * C[-2][k] - U[k] <= rows)  # Constraint 2e part 1
            if not isinstance(C[-1][k], (int, float)):
                model.addConstr((rows+1) * C[-1][k] + U[-1] - U[k] <= rows)  # Constraint 2e part 2
        
        # Add relation between dominating set and flow variables
        for i in range(rows):
            model.addConstr(var_1[i] == 1 - var_2[i])  # Constraint 2i
        
        # Set objective function
        z = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name='z')
        model.addConstr(z == sum(var_1))
        model.setObjective(z, cp.COPT.MINIMIZE)
        
        # Solve the model
        status = model.solve()
        
        # Process results
        if status == cp.COPT.OPTIMAL:
            VAR = [var.x for var in var_1]
            print(f'The connected domination number is: {model.objVal}')
        else:
            VAR = []
            print("No optimal solution found")
            
        return VAR
    
    def CMDS_SCF_OR(self):
        """Solve CMDS using OR-Tools CP-SAT solver with Single Commodity Flow formulation"""
        model = cp_model.CpModel()
        rows = self.nodes
        
        # Define variables
        var_1 = [model.NewBoolVar(name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_r = [model.NewBoolVar(name=f'r{i}') for i in range(rows)]  # Root variables
        var_e = [[model.NewIntVar(lb=0, ub=rows, name=f'f{i},{j}') for j in range(rows)] 
                for i in range(rows)]  # Flow variables
        
        # Compute matrix products
        B = self.adj * var_e
        B0 = self.B0 * var_1
        
        # Add domination constraints
        for i in range(rows):
            model.Add(sum(B0[i]) >= 1)
        
        # Define sum variables for easier expression
        sum_var_1 = model.NewIntVar(lb=0, ub=rows, name='sum_var_1')
        model.Add(sum_var_1 == cp_model.LinearExpr.Sum(var_1))
        
        # Add flow constraints
        var_z = []  # Multiplication results for var_r[i] * sum_var_1
        var_ml = []  # Multiplication results for var_1[i] * sum_var_1
        
        for i in range(rows):
            # Create variables for products
            z2 = model.NewIntVar(lb=0, ub=rows, name=f'z2_{i}')
            var_z.append(z2)
            model.AddMultiplicationEquality(z2, [var_r[i], sum_var_1])
            
            z3 = model.NewIntVar(lb=0, ub=rows, name=f'z3_{i}')
            var_ml.append(z3)
            model.AddMultiplicationEquality(z3, [var_1[i], sum_var_1])
            
            # Add membership constraints
            model.Add(var_1[i] >= var_r[i])
            
            # Add flow balance constraints
            sum_B_i = model.NewIntVar(lb=0, ub=rows, name=f'sum_B_{i}')
            model.Add(sum_B_i == cp_model.LinearExpr.Sum(B[i]))
            
            sum_B_T_i = model.NewIntVar(lb=0, ub=rows, name=f'sum_B_T_{i}')
            model.Add(sum_B_T_i == cp_model.LinearExpr.Sum([B[j][i] for j in range(rows)]))
            
            model.Add(sum_B_T_i <= rows * (1 - var_r[i]))
            model.Add(sum_B_T_i - sum_B_i - var_1[i] + var_z[i] == 0)
        
        # Add flow capacity constraints
        for i in range(rows):
            for j in range(rows):
                model.Add(B[i][j] <= var_ml[i])
                model.Add(B[i][j] <= var_ml[j])
        
        # Add root constraints
        model.Add(cp_model.LinearExpr.Sum(var_r) == 1)
        
        # Set objective function
        model.Minimize(cp_model.LinearExpr.Sum(var_1))
        
        # Solve the model
        solver = cp_model.CpSolver()
        status = solver.Solve(model)
        
        # Process results
        if status == cp_model.OPTIMAL:
            VAR = [solver.Value(var_1[i]) for i in range(rows)]
            print(f'The connected domination number is: {solver.ObjectiveValue()}')
        else:
            VAR = []
            print("No optimal solution found")
            
        return VAR
    
    def CMDS_SCF_GRB(self):
        """Solve CMDS using Gurobi solver with Single Commodity Flow formulation"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Define variables
        var_1 = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_r = [model.addVar(vtype=grb.GRB.BINARY, name=f'r{i}') for i in range(rows)]  # Root variables
        var_e = [[model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'f{i},{j}')
                for j in range(rows)] for i in range(rows)]  # Flow variables
        
        # Compute matrix products
        B = [[self.adj[j][k] * var_e[j][k] for k in range(rows)] for j in range(rows)]
        B0 = [[self.B0[j][k] * var_1[k] for k in range(rows)] for j in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(sum(B0[i]) >= 1)
        
        # Define sum variable
        sum_var_1 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name='sum_var_1')
        model.addConstr(sum_var_1 == sum(var_1))
        
        # Add flow constraints
        var_z = []  # Multiplication results for var_r[i] * sum_var_1
        var_ml = []  # Multiplication results for var_1[i] * sum_var_1
        
        for i in range(rows):
            # Create variables for products
            z2 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'z2_{i}')
            var_z.append(z2)
            model.addQConstr(z2 == var_r[i] * sum_var_1)
            
            z3 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'z3_{i}')
            var_ml.append(z3)
            model.addQConstr(z3 == var_1[i] * sum_var_1)
            
            # Add membership constraints
            model.addConstr(var_1[i] >= var_r[i])
            
            # Add flow balance constraints
            sum_B_i = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'sum_B_{i}')
            model.addConstr(sum_B_i == sum(B[i]))
            
            sum_B_T_i = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'sum_B_T_{i}')
            model.addConstr(sum_B_T_i == sum([B[j][i] for j in range(rows)]))
            
            model.addConstr(sum_B_T_i <= rows * (1 - var_r[i]))
            model.addConstr(sum_B_T_i - sum_B_i - var_1[i] + var_z[i] == 0)
        
        # Add flow capacity constraints
        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j] <= var_ml[i])
                model.addConstr(B[i][j] <= var_ml[j])
        
        # Add root constraints
        model.addConstr(sum(var_r) == 1)
        
        # Set objective function
        model.setObjective(sum(var_1), grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        VAR = [var_1[i].x for i in range(rows)]
        print(f'The connected domination number is: {model.objVal}')
        
        return VAR
    
    def CMDS_MARTIN_GRB(self):
        """Solve CMDS using Gurobi solver with Martin formulation"""
        model = self._create_gurobi_env()
        model.setParam('TimeLimit', 5*60)  # 5 minutes time limit
        rows = self.nodes
        M = self.M
        
        # Define variables
        var_x = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_y = [[model.addVar(vtype=grb.GRB.BINARY, name=f'y{i}{j}') for j in range(rows)] 
                for i in range(rows)]  # Edge variables
        var_z = [[[model.addVar(vtype=grb.GRB.BINARY, name=f'z{i}{j}{k}') for k in range(rows)]
                for j in range(rows)] for i in range(rows)]  # Path variables
        
        # Compute matrix products
        A0 = [[self.adj[i][j] * var_x[j] for j in range(rows)] for i in range(rows)]
        B = [[self.adj[i][j] * var_y[i][j] for j in range(rows)] for i in range(rows)]
        C = [[[self.adj[i][j] * var_z[i][j][k] for k in range(rows)] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(sum(A0[i]) >= 1)
        
        # Add connectivity constraints
        model.addConstr(sum([sum(B[i]) for i in range(rows)]) == sum(var_x) - 1)
        
        # Add flow constraints
        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j] <= var_x[i])
                model.addConstr(B[i][j] <= var_x[j])
                
                for k in range(rows):
                    model.addConstr(C[i][j][k] <= B[i][j])
                    model.addConstr(C[i][j][k] <= var_x[k])
                    model.addConstr(C[j][i][k] <= B[i][j])
                    model.addConstr(C[j][i][k] <= var_x[k])
                    model.addConstr(B[i][j] - M * (3 - var_x[i] - var_x[j] - var_x[k]) - C[i][j][k] - C[j][i][k] <= 0)
                    model.addConstr(C[i][j][k] + C[j][i][k] <= B[i][j] + M * (3 - var_x[i] - var_x[j] - var_x[k]))
        
        # Add path constraints
        D = [[[C[i][j][k] for j in range(rows) if j != k and j != i] for k in range(rows)] for i in range(rows)]
        for i in range(rows):
            for k in range(rows):
                model.addConstr(sum(D[i][k]) + B[i][k] >= 1 - M * (2 - var_x[i] - var_x[j]))
                model.addConstr(sum(D[i][k]) + B[i][k] <= 1 + M * (2 - var_x[i] - var_x[j]))
        
        # Set objective function
        model.setObjective(sum(var_x), grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        VAR = [var_x[i].x for i in range(rows)]
        print(f'The connected domination number is: {model.objVal}')
        
        return VAR
    
    def CMDS_MARTIN_CP(self):
        """Solve CMDS using COPT solver with Martin formulation"""
        model = self._create_copt_env()
        model.setParam(cp.COPT.Param.TimeLimit, 5*60)  # 5 minutes time limit
        rows = self.nodes
        M = self.M
        
        # Define variables
        var_x = [model.addVar(vtype=cp.COPT.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_y = [[model.addVar(vtype=cp.COPT.BINARY, name=f'y{i}{j}') for j in range(rows)] 
                for i in range(rows)]  # Edge variables
        var_z = [[[model.addVar(vtype=cp.COPT.BINARY, name=f'z{i}{j}{k}') for k in range(rows)]
                for j in range(rows)] for i in range(rows)]  # Path variables
        
        # Compute matrix products
        A0 = [[self.adj[i][j] * var_x[j] for j in range(rows)] for i in range(rows)]
        B = [[self.adj[i][j] * var_y[i][j] for j in range(rows)] for i in range(rows)]
        C = [[[self.adj[i][j] * var_z[i][j][k] for k in range(rows)] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(sum(A0[i]) >= 1)
        
        # Add connectivity constraints
        model.addConstr(sum([sum(B[i]) for i in range(rows)]) == sum(var_x) - 1)
        
        # Add flow constraints
        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j] <= var_x[i])
                model.addConstr(B[i][j] <= var_x[j])
                
                for k in range(rows):
                    model.addConstr(C[i][j][k] <= B[i][j])
                    model.addConstr(C[i][j][k] <= var_x[k])
                    model.addConstr(C[j][i][k] <= B[i][j])
                    model.addConstr(C[j][i][k] <= var_x[k])
                    model.addConstr(B[i][j] - M * (3 - var_x[i] - var_x[j] - var_x[k]) - C[i][j][k] - C[j][i][k] <= 0)
                    model.addConstr(C[i][j][k] + C[j][i][k] <= B[i][j] + M * (3 - var_x[i] - var_x[j] - var_x[k]))
        
        # Add path constraints
        D = [[[C[i][j][k] for j in range(rows) if j != k and j != i] for k in range(rows)] for i in range(rows)]
        for i in range(rows):
            for k in range(rows):
                model.addConstr(sum(D[i][k]) + B[i][k] >= 1 - M * (2 - var_x[i] - var_x[j]))
                model.addConstr(sum(D[i][k]) + B[i][k] <= 1 + M * (2 - var_x[i] - var_x[j]))
        
        # Set objective function
        model.setObjective(sum(var_x), cp.COPT.MINIMIZE)
        
        # Solve the model
        stat = model.solve()
        
        # Process results
        if stat == cp.COPT.OPTIMAL:
            VAR = [var_x[i].X for i in range(rows)]
            print(f'The connected domination number is: {model.objVal}')
        else:
            VAR = []
            print('No feasible solution.')
        
        return VAR
    
    def CMDS(self, cons):
        """Solve CMDS using the selected solver and constraint formulation"""
        solver = self.solver
        if solver == 'OR':
            if cons in ['MTZ', 'SCF']:
                return getattr(self, f'CMDS_{cons}_OR')()
            else:
                raise ValueError(f"Constraint formulation {cons} not supported for OR-Tools")
        elif solver == 'GRB':
            if cons in ['MTZ', 'SCF', 'MAR']:
                return getattr(self, f'CMDS_{cons}_GRB')()
            else:
                raise ValueError(f"Constraint formulation {cons} not supported for Gurobi")
        elif solver == 'CP':
            if cons == 'SCF':
                print('Sorry! COPT cannot solve this problem')
                return []
            elif cons in ['MTZ', 'MAR']:
                return getattr(self, f'CMDS_{cons}_CP')()
            else:
                raise ValueError(f"Constraint formulation {cons} not supported for COPT")
        else:
            raise ValueError(f"Unknown solver: {self.solver}")
    
    #-------------------------------------------------------------------------
    # Power Minimum Dominating Set (PMDS) Methods
    #-------------------------------------------------------------------------
    
    def PMDS_GRB(self, Z):
        """Solve PMDS using Gurobi solver"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Define variables
        var_1 = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_p = [[model.addVar(vtype=grb.GRB.BINARY, name=f'p{i}{j}') for j in range(rows)] 
                for i in range(rows)]  # Power variables
        
        # Add power constraints
        for i in range(rows):
            for j in range(rows):
                if j != i:
                    model.addConstr(var_p[i][j] <= self.adj[i][j])
                else:
                    model.addConstr(var_p[i][j] == 0)
        
        # Compute matrix products
        B = [[self.adj[i][j] * var_1[j] for j in range(rows)] for i in range(rows)]
        C = [[self.adj[i][j] * var_p[i][j] for j in range(rows)] for i in range(rows)]
        Z0 = [[self.adj[i][j] * Z[j] for j in range(rows)] for i in range(rows)]
        Z1 = [[Z0[i][j] * var_p[j][i] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(sum(B[i]) + sum(Z1[i]) >= 1)
            model.addConstr(sum(C[i]) == Z[i])
        
        # Set objective function
        model.setObjective(sum(var_1), grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        VAR = [var_1[i].x for i in range(rows)]
        print(f'The power domination number is: {model.objVal}')
        
        return VAR
    
    def PMDS_CP(self, Z):
        """Solve PMDS using COPT solver"""
        model = self._create_copt_env()
        rows = self.nodes
        
        # Define variables
        var_1 = [model.addVar(vtype=cp.COPT.BINARY, name=f'x{i}') for i in range(rows)]  # Dominating set variables
        var_p = [[model.addVar(vtype=cp.COPT.BINARY, name=f'p{i}{j}') for j in range(rows)] 
                for i in range(rows)]  # Power variables
        
        # Add power constraints
        for i in range(rows):
            for j in range(rows):
                if j != i:
                    model.addConstr(var_p[i][j] <= self.adj[i][j])
                else:
                    model.addConstr(var_p[i][j] == 0)
        
        # Compute matrix products
        B = [[self.adj[i][j] * var_1[j] for j in range(rows)] for i in range(rows)]
        C = [[self.adj[i][j] * var_p[i][j] for j in range(rows)] for i in range(rows)]
        Z0 = [[self.adj[i][j] * Z[j] for j in range(rows)] for i in range(rows)]
        Z1 = [[Z0[i][j] * var_p[j][i] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(sum(B[i]) + sum(Z1[i]) >= 1)
            model.addConstr(sum(C[i]) == Z[i])
        
        # Set objective function
        model.setObjective(sum(var_1), cp.COPT.MINIMIZE)
        
        # Solve the model
        status = model.solve()
        
        # Process results
        if status == cp.COPT.OPTIMAL:
            VAR = [var_1[i].x for i in range(rows)]
            print(f'The power domination number is: {model.objVal}')
        else:
            VAR = []
            print("No optimal solution found")
        
        return VAR
    
    def PMDS(self, Z):
        """Solve PMDS using the selected solver"""
        if self.solver == 'GRB':
            return self.PMDS_GRB(Z)
        elif self.solver == 'CP':
            return self.PMDS_CP(Z)
        else:
            raise ValueError(f"Solver {self.solver} not supported for PMDS")
    
    #-------------------------------------------------------------------------
    # Distance-Related Domination Problem (DRDP) Methods
    #-------------------------------------------------------------------------
    
    def _solve_drdp(self, model_type, solver):
        """Helper method to solve DRDP with the specified model and solver"""
        method_name = f'DRDP_{model_type}_{solver}'
        if hasattr(self, method_name):
            return getattr(self, method_name)()
        else:
            raise ValueError(f"Model {model_type} not supported for solver {solver}")
    
    def DRDP_1_GRB(self):
        """Solve DRDP model 1 using Gurobi solver"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Define variables
        var_x = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]  # Distance-1 variables
        var_y = [model.addVar(vtype=grb.GRB.BINARY, name=f'y{i}') for i in range(rows)]  # Distance-2 variables
        var_z = [model.addVar(vtype=grb.GRB.BINARY, name=f'z{i}') for i in range(rows)]  # Distance-3 variables
        
        # Compute matrix products
        B = [[self.adj[i][j] * var_y[j] for j in range(rows)] for i in range(rows)]
        C = [[self.adj[i][j] * var_z[j] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(var_x[i] + var_y[i] + var_z[i] + 0.5 * sum(B[i]) + sum(C[i]) >= 1)
            model.addConstr(var_x[i] <= sum(B[i]) + sum(C[i]))
            model.addConstr(var_x[i] + var_y[i] + var_z[i] <= 1)
        
        # Set objective function
        model.setObjective(sum(var_x) + 2 * sum(var_y) + 3 * sum(var_z), grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        ans = []
        for i in range(rows):
            if var_x[i].x == 1:
                ans.append(1)
            elif var_y[i].x == 1:
                ans.append(2)
            elif var_z[i].x == 1:
                ans.append(3)
            else:
                ans.append(0)
        
        print(f'DRDP model 1 result: {ans}')
        return ans
    
    def DRDP_1_CP(self):
        """Solve DRDP model 1 using COPT solver"""
        model = self._create_copt_env()
        rows = self.nodes
        
        # Define variables
        var_x = [model.addVar(vtype=cp.COPT.BINARY, name=f'x{i}') for i in range(rows)]  # Distance-1 variables
        var_y = [model.addVar(vtype=cp.COPT.BINARY, name=f'y{i}') for i in range(rows)]  # Distance-2 variables
        var_z = [model.addVar(vtype=cp.COPT.BINARY, name=f'z{i}') for i in range(rows)]  # Distance-3 variables
        
        # Compute matrix products
        B = [[self.adj[i][j] * var_y[j] for j in range(rows)] for i in range(rows)]
        C = [[self.adj[i][j] * var_z[j] for j in range(rows)] for i in range(rows)]
        
        # Add domination constraints
        for i in range(rows):
            model.addConstr(var_x[i] + var_y[i] + var_z[i] + 0.5 * sum(B[i]) + sum(C[i]) >= 1)
            model.addConstr(var_x[i] <= sum(B[i]) + sum(C[i]))
            model.addConstr(var_x[i] + var_y[i] + var_z[i] <= 1)
        
        # Set objective function
        model.setObjective(sum(var_x) + 2 * sum(var_y) + 3 * sum(var_z), cp.COPT.MINIMIZE)
        
        # Solve the model
        model.solve()
        
        # Process results
        ans = []
        for i in range(rows):
            if var_x[i].x == 1:
                ans.append(1)
            elif var_y[i].x == 1:
                ans.append(2)
            elif var_z[i].x == 1:
                ans.append(3)
            else:
                ans.append(0)
        
        print(f'DRDP model 1 result: {ans}')
        return ans
    
    def DRDP_2_GRB(self):
        """Solve DRDP model 2 using Gurobi solver"""
        model = self._create_gurobi_env()
        rows = self.nodes
        
        # Define variables
        var_p = [model.addVar(vtype=grb.GRB.BINARY, name=f'p{i}') for i in range(rows)]
        var_q = [model.addVar(vtype=grb.GRB.BINARY, name=f'q{i}') for i in range(rows)]
        var_r = [model.addVar(vtype=grb.GRB.BINARY, name=f'r{i}') for i in range(rows)]
        
        # Compute matrix products
        B = [[self.adj[i][j] * var_q[j] for j in range(rows)] for i in range(rows)]
        C = [[self.adj[i][j] * var_r[j] for j in range(rows)] for i in range(rows)]
        
        # Add constraints
        for i in range(rows):
            model.addConstr(var_p[i] + 0.5 * sum(B[i]) + 0.5 * sum(C[i]) >= 1)
            model.addConstr(var_q[i] + sum(B[i]) >= var_p[i])
            model.addConstr(var_p[i] >= var_q[i])
            model.addConstr(var_q[i] >= var_r[i])
        
        # Set objective function
        model.setObjective(sum(var_p) + sum(var_q) + sum(var_r), grb.GRB.MINIMIZE)
        
        # Solve the model
        model.optimize()
        
        # Process results
        ans = []
        for i in range(rows):
            if var_p[i].x == 1 and var_q[i].x == 0:
                ans.append(1)
            elif var_q[i].x == 1 and var_p[i].x == 1:
                ans.append(2)
            elif var_r[i].x == 1:
                ans.append(3)
            else:
                ans.append(0)
        
        print(f'DRDP model 2 result: {ans}')
        return ans
    
    def DRDP(self, n):
        """Solve DRDP using the selected solver and model"""
        if self.solver == 'GRB':
            return self._solve_drdp(n, 'GRB')
        elif self.solver == 'CP':
            return self._solve_drdp(n, 'CP')
        else:
            raise ValueError(f"Solver {self.solver} not supported for DRDP")
    
    #-------------------------------------------------------------------------
    # Visualization Methods
    #-------------------------------------------------------------------------
    
    def draw(self, var):
        """Visualize the graph with highlighted dominating nodes"""
        # Construct the graph using the adjacency matrix
        G = nx.from_numpy_array(self.adj)
        
        # Set labels for the nodes
        node_labels = {i: str(i) for i in range(self.nodes)}
        
        # Determine which nodes are in the dominating set
        dom = [i for i in range(self.nodes) if var[i] == 1]
        
        # Draw the graph
        plt.figure(figsize=(16, 8))
        colors = ['r' if i in dom else 'b' for i in range(self.nodes)]
        nx.draw(G, labels=node_labels, with_labels=True, node_color=colors)
        plt.title("Dominating Set Visualization")
        plt.show()
    
    def draw_DRDP(self, var):
        """Visualize the graph with distance-related domination coloring"""
        # Construct the graph using the adjacency matrix
        G = nx.from_numpy_array(self.adj)
        
        # Set labels for the nodes
        node_labels = {i: str(i) for i in range(self.nodes)}
        
        # Define colors for different distance values
        colors = ['r', 'b', 'y', 'g']
        col = [colors[var[i]] for i in range(self.nodes)]
        
        # Draw the graph
        plt.figure(figsize=(16, 8))
        nx.draw(G, labels=node_labels, with_labels=True, node_color=col)
        plt.title("Distance-Related Domination Visualization")
        plt.show()


# Rename the class for backward compatibility
dom = DominationProblemSolver