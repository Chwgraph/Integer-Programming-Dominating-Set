# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:10:05 2025

@author: howiech
"""

import numpy as np
import networkx as nx
import gurobipy as grb
from ortools.sat.python import cp_model
import coptpy as cp
from coptpy import COPT
import matplotlib.pyplot as plt



class dom:
    def __init__(self, adj_matrix, nodes, solver, params, M):
        self.adj=adj_matrix
        self.nodes=nodes
        self.solver=solver
        self.params=params
        self.M=M
    
    #Solver MDS problem using CP solver in OR_TOOLS
    def MDS_OR(self, size, version):
        model=cp_model.CpModel()
        A = np.array(self.adj)
        rows=self.nodes
        B=A+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        #Add variables to the model
        variables = []
        for i in range(rows):
            var = model.NewIntVar(0, 1, 'x'+str(i))
            variables.append(var)
        #Add target variables
        z = model.NewIntVar(0, rows, 'z')

        # For each node, either at least one of its neighbors is in the dominating set
        C = B*variables
        for j in range(rows):
            if version =='w':
                if C[j][j]==0:
                    model.Add(sum(C[j])>=size)
            elif version == 's':
                model.Add(sum(C[j])>=size)
            
        model.Add(z==sum(variables))

        #Set objective function
        model.Minimize(z)
        #Initialize the solver
        solver=cp_model.CpSolver()

        #Solve the model
        status=solver.Solve(model)

        #If the solver finds the optimal solution, print the solution and the value of the objective function
        if status==cp_model.OPTIMAL:
            for var in variables:
                print(solver.Value(var))
            print('The domination number of the input graph is:{}'.format(solver.Value(z)))
            print('The number of vertices in the input graph is{}'.format(rows))
        VAR=[solver.Value(var) for var in variables]
        return VAR
    
    #Solve MDS using Gurobi Solver
    def MDS_GRB(self, size, version):
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        A = np.array(self.adj)
        rows=self.nodes
        B=A+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        

        #Add variables to the model
        variables=[]
        for i in range(rows):
            var=model.addVar(vtype=grb.GRB.BINARY,lb=0,ub=2,name='x'+str(i))
            variables.append(var)

        #Add the objective variable
        z=model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name='z')

        # For each node, either at least one of its neighbors is in the dominating set
        C=B*variables

        for j in range(rows):
            if version == 'w':
                if C[j][j]==0:
                    model.addConstr(sum(C[j])>=size)
            if version == 's':
                model.addConstr(sum(C[j])>=size)
            
        #Set the objective function
        model.addConstr(z==sum(variables))
        model.setObjective(z,grb.GRB.MINIMIZE)

        #Solve the model
        model.optimize()


        VAR=[var.x for var in variables]
        print('The {}-domination number is:{}'.format(size, model.objval))
        print(VAR)
        return VAR
    
    #Solve MDS using COPT Solver
    def MDS_CP(self, size, version):
        #Import the solver environment
        env=cp.Envr()

        #Create the model
        model=env.createModel('DOM')
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        A = np.array(self.adj)
        rows=self.nodes
        B=A+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])

        #Add the variables to the model
        variables=[]
        for i in range(rows):
            var=model.addVar(vtype=cp.COPT.BINARY, lb=0, ub=2, name='x'+str(i))
            variables.append(var)

        #Add the objective variable
        z=model.addVar(vtype=cp.COPT.INTEGER, name='z')

        # For each node, either at least one of its neighbors is in the dominating set
        C=[[B[j][k]*variables[k] for k in range(rows)] for j in range(rows)]
            
        for j in range(rows):
            if version == 'w':
                if C[j][j] == 0:
                    model.addConstr(sum(C[j])>=size)
            if version == 's':
                model.addConstr(sum(C[j])>=size)

        #Set the objective function
        model.addConstr(z==sum(variables))
        model.setObjective(z, cp.COPT.MINIMIZE)

        #Solve the model
        status=model.solve()

        #If the solver finds the optimal solution, print the solution
        print('The {}-domination number is:{}'.format(size, model.objVal))
        VAR=[var.x for var in variables]
        print(VAR)
        return VAR
    
    #Solver MDS using the solver selectedby the user.
    def MDS(self, k, ver):
        solver=self.solver
        if solver=='OR':
            var=self.MDS_OR(k, ver)
        elif solver=='GRB':
            var=self.MDS_GRB(k, ver)
        elif solver=='CP':
            var=self.MDS_CP(k, ver)
        return var
    
    def CMDS_MTZ_OR(self):
        model=cp_model.CpModel()
        rows=self.nodes
        A = [list(self.adj[i]) for i in range(rows)]+[[1 for i in range(rows)], [1 for i in range(rows)]]
        B=np.array(A)
        B0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        var_1=[]
        var_2=[]
        var_3=[]
        var_4=[]
        U=[]
        for i in range(rows):
            var1=model.NewBoolVar(name='x'+str(i))
            var2=model.NewBoolVar(name=r'yn+1,{}'.format(i))
            var3=model.NewBoolVar(name=r'yn+2,{}'.format(i))
            u=model.NewIntVar(1,rows+1, name=r'u{}'.format(i))
            var_1.append(var1)
            var_2.append(var2)
            var_3.append(var3)
            U.append(u)
            var_row=[]
            for j in range(rows):
                var4=model.NewBoolVar(name=r'y{},{}'.format(i,j))
                var_row.append(var4)
            var_4.append(var_row)
        var_4.append(var_2)
        var_4.append(var_3)
        U.append(0)
        u=model.NewIntVar(1, rows+1, name=r'u{}'.format(rows+2))
        U.append(u)
                    
        C=np.array(B*var_4)
        C0=B0*var_1
        for j in range(rows):
            model.Add(sum(C0[j])>=1)
        #2a
        model.Add(sum(var_3)==1)
        #2b
        for k in range(rows):
            model.Add(sum(C.T[k])==1)
        #2c
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.Add(var_2[k]+C[k][l]<=1)
        
        #2d
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.Add((rows+1)*C[k][l]+U[k]-U[l]+(rows-1)*C[l][k]<=rows)
        
        #2e
        for k in range(rows):
            if type(C[-2][k])!=int:
                model.Add((rows+1)*C[-2][k]-U[k]<=rows)
            if type(C[-1][k])!=int:
                model.Add((rows+1)*C[-1][k]+U[-1]-U[k]<=rows)
        
        #2i
        for i in range(rows):
            model.Add(var_1[i]==1-var_2[i])
        
        z=model.NewIntVar(lb=0,ub=rows,name='z')
        model.Add(z==sum(var_1))
        
        #Set the objective function.
        model.Minimize(z)
        
        #Solve the model
        solver=cp_model.CpSolver()
        status=solver.Solve(model)
        if status==cp_model.OPTIMAL:
            VAR=[solver.Value(var_1[i]) for i in range(rows)]
        print(VAR)
        return VAR
    
    def CMDS_MTZ_GRB(self):
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        rows=self.nodes
        A = [list(self.adj[i]) for i in range(rows)]+[[1 for i in range(rows)], [1 for i in range(rows)]]
        B=np.array(A)
        B0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        var_1=[]
        var_2=[]
        var_3=[]
        var_4=[]
        U=[]
        for i in range(rows):
            var1=model.addVar(vtype=grb.GRB.BINARY,name='x'+str(i))
            var2=model.addVar(vtype=grb.GRB.BINARY, name=r'yn+1,{}'.format(i))
            var3=model.addVar(vtype=grb.GRB.BINARY, name=r'yn+2,{}'.format(i))
            u=model.addVar(vtype=grb.GRB.INTEGER, lb=1,ub=rows+1, name=r'u{}'.format(i))
            var_1.append(var1)
            var_2.append(var2)
            var_3.append(var3)
            U.append(u)
            var_row=[]
            for j in range(rows):
                var4=model.addVar(vtype=grb.GRB.BINARY, name=r'y{},{}'.format(i,j))
                var_row.append(var4)
            var_4.append(var_row)
        var_4.append(var_2)
        var_4.append(var_3)
        U.append(0)
        u=model.addVar(vtype=grb.GRB.INTEGER, lb=1, ub=rows+1, name=r'u{}'.format(rows+2))
        U.append(u)
                    
        C=np.array(B*var_4)
        C0=B0*var_1
        for j in range(rows):
            model.addConstr(sum(C0[j])>=1)
        #2a
        model.addConstr(sum(var_3)==1)
        #2b
        for k in range(rows):
            model.addConstr(sum(C.T[k])==1)
        #2c
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.addConstr(var_2[k]+C[k][l]<=1)
        
        #2d
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.addConstr((rows+1)*C[k][l]+U[k]-U[l]+(rows-1)*C[l][k]<=rows)
        
        #2e
        for k in range(rows):
            if type(C[-2][k])!=int:
                model.addConstr((rows+1)*C[-2][k]-U[k]<=rows)
            if type(C[-1][k])!=int:
                model.addConstr((rows+1)*C[-1][k]+U[-1]-U[k]<=rows)
        
        #2i
        for i in range(rows):
            model.addConstr(var_1[i]==1-var_2[i])
        
        z=model.addVar(vtype=grb.GRB.INTEGER, lb=0,ub=rows,name='z')
        model.addConstr(z==sum(var_1))
        
        #Set the objective function.
        model.setObjective(z, grb.GRB.MINIMIZE)
        model.optimize()
        print('The connected domination number is:{}'.format(model.objVal))
        VAR=[var.x for var in var_1]
        print(VAR)
        return VAR
    
    def CMDS_MTZ_CP(self):
        env=cp.Envr()
        model=env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        rows=self.nodes
        A = [list(self.adj[i]) for i in range(rows)]+[[1 for i in range(rows)], [1 for i in range(rows)]]
        B=np.array(A)
        B0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        var_1=[]
        var_2=[]
        var_3=[]
        var_4=[]
        U=[]
        for i in range(rows):
            var1=model.addVar(vtype=cp.COPT.BINARY,name='x'+str(i))
            var2=model.addVar(vtype=cp.COPT.BINARY, name=r'yn+1,{}'.format(i))
            var3=model.addVar(vtype=cp.COPT.BINARY, name=r'yn+2,{}'.format(i))
            u=model.addVar(vtype=cp.COPT.INTEGER, lb=1,ub=rows+1, name=r'u{}'.format(i))
            var_1.append(var1)
            var_2.append(var2)
            var_3.append(var3)
            U.append(u)
            var_row=[]
            for j in range(rows):
                var4=model.addVar(vtype=cp.COPT.BINARY, name=r'y{},{}'.format(i,j))
                var_row.append(var4)
            var_4.append(var_row)
        var_4.append(var_2)
        var_4.append(var_3)
        U.append(0)
        u=model.addVar(vtype=cp.COPT.INTEGER, lb=1, ub=rows+1, name=r'u{}'.format(rows+2))
        U.append(u)
                    
        C=np.array([[B[j][k]*var_4[j][k] for k in range(rows)] for j in range(rows+2)])
        print(C)
        C0=[[B0[j][k]*var_1[k] for k in range(rows)] for j in range(rows)]
        for j in range(rows):
            model.addConstr(sum(C0[j])>=1)
        #2a
        model.addConstr(sum(var_3)==1)
        #2b
        for k in range(rows):
            model.addConstr(sum(C.T[k])==1)
        #2c
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.addConstr(var_2[k]+C[k][l]<=1)
        
        #2d
        for k in range(rows):
            for l in range(rows):
                if type(C[k][l])!=int:
                    model.addConstr((rows+1)*C[k][l]+U[k]-U[l]+(rows-1)*C[l][k]<=rows)
        
        #2e
        for k in range(rows):
            if type(C[-2][k])!=int:
                model.addConstr((rows+1)*C[-2][k]-U[k]<=rows)
            if type(C[-1][k])!=int:
                model.addConstr((rows+1)*C[-1][k]+U[-1]-U[k]<=rows)
        
        #2i
        for i in range(rows):
            model.addConstr(var_1[i]==1-var_2[i])
        
        z=model.addVar(vtype=cp.COPT.INTEGER, lb=0,ub=rows,name='z')
        model.addConstr(z==sum(var_1))
        
        #Set the objective function.
        model.setObjective(z, cp.COPT.MINIMIZE)
        model.solve()
        print('The connected domination number is:{}'.format(model.objVal))
        VAR=[var.x for var in var_1]
        return VAR
    
    def CMDS_SCF_OR(self):
        model = cp_model.CpModel()
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])

        # Define variables
        var_1 = [model.NewBoolVar(name=f'x{i}') for i in range(rows)]
        var_r = [model.NewBoolVar(name=f'r{i}') for i in range(rows)]
        var_e = [[model.NewIntVar(lb=0, ub=rows, name=f'f{i},{j}') for j in range(rows)] for i in range(rows)]

        # Compute B = A * var_e
        B=A*var_e
        B0=A0*var_1
        for i in range(rows):
            model.Add(sum(B0[i])>=1)

        # Define sum variables
        sum_var_1 = model.NewIntVar(lb=0, ub=rows, name='sum_var_1')
        model.Add(sum_var_1 == cp_model.LinearExpr.Sum(var_1))

        var_z = []
        var_ml = []
        for i in range(rows):
            z2 = model.NewIntVar(lb=0, ub=rows, name=f'z2_{i}')
            var_z.append(z2)
            model.AddMultiplicationEquality(z2, [var_r[i], sum_var_1])

            z3 = model.NewIntVar(lb=0, ub=rows, name=f'z3_{i}')
            var_ml.append(z3)
            model.AddMultiplicationEquality(z3, [var_1[i], sum_var_1])

            model.Add(var_1[i] >= var_r[i])

            sum_B_i = model.NewIntVar(lb=0, ub=rows, name=f'sum_B_{i}')
            model.Add(sum_B_i == cp_model.LinearExpr.Sum(B[i]))

            sum_B_T_i = model.NewIntVar(lb=0, ub=rows, name=f'sum_B_T_{i}')
            model.Add(sum_B_T_i == cp_model.LinearExpr.Sum([B[j][i] for j in range(rows)]))

            model.Add(sum_B_T_i <= rows * (1 - var_r[i]))
            model.Add(sum_B_T_i - sum_B_i - var_1[i] + var_z[i] == 0)

        for i in range(rows):
            for j in range(rows):
                model.Add(B[i][j] <= var_ml[i])
                model.Add(B[i][j] <= var_ml[j])

        model.Add(cp_model.LinearExpr.Sum(var_r) == 1)
        model.Minimize(cp_model.LinearExpr.Sum(var_1))

        # Solve the model
        solver = cp_model.CpSolver()
        status = solver.Solve(model)

        # Extract results
        VAR = [solver.Value(var_1[i]) for i in range(rows)]
        if status == cp_model.OPTIMAL:
            print(VAR)
        return VAR
    
    # def CMDS_SCF_CP(self):
    #     rows = self.nodes
    #     A = np.array(self.adj)
    #     A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
    #     envr=cp.Envr()
    #     model=envr.createModel()
        
    #     # Define variables
    #     var_1 = [model.addVar(vtype=cp.COPT.BINARY, name=f'x{i}') for i in range(rows)]
    #     var_r = [model.addVar(vtype=cp.COPT.BINARY, name=f'r{i}') for i in range(rows)]
    #     var_e = [[model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name=f'f{i},{j}') for j in range(rows)] for i in range(rows)]

    #     # Compute B = A * var_e
    #     B=[[A[j][k]*var_e[j][k] for k in range(rows)] for j in range(rows)]
    #     B0=[[A0[j][k]*var_1[k] for k in range(rows)] for j in range(rows)]
    #     for i in range(rows):
    #         model.addConstr(sum(B0[i])>=1)

    #     # Define sum variables
    #     sum_var_1 = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name='sum_var_1')
    #     model.addConstr(sum_var_1 == sum(var_1))

    #     var_z = []
    #     var_ml = []
    #     for i in range(rows):
    #         z2 = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name=f'z2_{i}')
    #         var_z.append(z2)
    #         model.addQConstr(z2==var_r[i]*sum_var_1)

    #         z3 = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name=f'z3_{i}')
    #         var_ml.append(z3)
    #         model.addQConstr(z3==var_1[i]*sum_var_1)

    #         model.addConstr(var_1[i] >= var_r[i])

    #         sum_B_i = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name=f'sum_B_{i}')
    #         model.addConstr(sum_B_i == sum(B[i]))

    #         sum_B_T_i = model.addVar(vtype=cp.COPT.INTEGER, lb=0, ub=rows, name=f'sum_B_T_{i}')
    #         model.addConstr(sum_B_T_i == sum([B[j][i] for j in range(rows)]))

    #         model.addConstr(sum_B_T_i <= rows * (1 - var_r[i]))
    #         model.addConstr(sum_B_T_i - sum_B_i - var_1[i] + var_z[i] == 0)

    #     for i in range(rows):
    #         for j in range(rows):
    #             model.addConstr(B[i][j] <= var_ml[i])
    #             model.addConstr(B[i][j] <= var_ml[j])

    #     model.addConstr(sum(var_r) == 1)
    #     model.setObjective(sum(var_1), cp.COPT.MINIMIZE)

    #     # Solve the model
    #     status = model.solve()
    #     # Extract results
    #     VAR = [var_1[i].x for i in range(rows)]
    #     if status == cp.COPT.OPTIMAL:
    #         print(VAR)
    #     return VAR
    
    def CMDS_SCF_GRB(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        
        # Define variables
        var_1 = [model.addVar(vtype=grb.GRB.BINARY, name=f'x{i}') for i in range(rows)]
        var_r = [model.addVar(vtype=grb.GRB.BINARY, name=f'r{i}') for i in range(rows)]
        var_e = [[model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'f{i},{j}') for j in range(rows)] for i in range(rows)]

        # Compute B = A * var_e
        B=[[A[j][k]*var_e[j][k] for k in range(rows)] for j in range(rows)]
        B0=[[A0[j][k]*var_1[k] for k in range(rows)] for j in range(rows)]
        for i in range(rows):
            model.addConstr(sum(B0[i])>=1)

        # Define sum variables
        sum_var_1 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name='sum_var_1')
        model.addConstr(sum_var_1 == sum(var_1))

        var_z = []
        var_ml = []
        for i in range(rows):
            z2 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'z2_{i}')
            var_z.append(z2)
            model.addQConstr(z2==var_r[i]*sum_var_1)

            z3 = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'z3_{i}')
            var_ml.append(z3)
            model.addQConstr(z3==var_1[i]*sum_var_1)

            model.addConstr(var_1[i] >= var_r[i])

            sum_B_i = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'sum_B_{i}')
            model.addConstr(sum_B_i == sum(B[i]))

            sum_B_T_i = model.addVar(vtype=grb.GRB.INTEGER, lb=0, ub=rows, name=f'sum_B_T_{i}')
            model.addConstr(sum_B_T_i == sum([B[j][i] for j in range(rows)]))

            model.addConstr(sum_B_T_i <= rows * (1 - var_r[i]))
            model.addConstr(sum_B_T_i - sum_B_i - var_1[i] + var_z[i] == 0)

        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j] <= var_ml[i])
                model.addConstr(B[i][j] <= var_ml[j])

        model.addConstr(sum(var_r) == 1)
        model.setObjective(sum(var_1), grb.GRB.MINIMIZE)

        # Solve the model
        status = model.optimize()
        # Extract results
        VAR = [var_1[i].x for i in range(rows)]
        if status == grb.GRB.OPTIMAL:
            print(VAR)
        return VAR
        
    def CMDS_MARTIN_GRB(self):
        rows = self.nodes
        M=self.M
        A = np.array(self.adj)
        A0 = np.array(self.adj) + np.array([[0 if j != i else 1 for j in range(rows)] for i in range(rows)])
        params = self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 5 * 60)

        var_x=[]
        var_y=[]
        var_z=[]
        for i in range(rows):
            x=model.addVar(vtype=grb.GRB.BINARY, name='x'+str(i))
            var_x.append(x)
            var_1=[]
            var_2=[]
            for j in range(rows):
                y=model.addVar(vtype=grb.GRB.BINARY, name='y'+str(i)+str(j))
                var_1.append(y)
                var_3=[]
                for k in range(rows):
                    z=model.addVar(vtype=grb.GRB.BINARY, name='z'+str(i)+str(j)+str(k))
                    var_3.append(z)
                var_2.append(var_3)
            var_y.append(var_1)
            var_z.append(var_2)
        A0 = [[A[i][j] * var_x[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(sum(A0[i])>=1)
        B=[[A[i][j]*var_y[i][j] for j in range(rows)] for i in range(rows)]
        C=[[[A[i][j]*var_z[i][j][k] for k in range(rows)] for j in range(rows)] for i in range(rows)]
        model.addConstr(sum([sum(B[i]) for i in range(rows)])==sum(var_x)-1)
        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j]<=var_x[i])
                model.addConstr(B[i][j]<=var_x[j])
                for k in range(rows):
                    model.addConstr(C[i][j][k]<=B[i][j])
                    model.addConstr(C[i][j][k]<=var_x[k])
                    model.addConstr(C[j][i][k]<=B[i][j])
                    model.addConstr(C[j][i][k] <= var_x[k])
                    model.addConstr(B[i][j]-M*(3-var_x[i]-var_x[j]-var_x[k])-C[i][j][k]-C[j][i][k]<=0)
                    model.addConstr(C[i][j][k]+C[j][i][k]<=B[i][j]+M*(3-var_x[i]-var_x[j]-var_x[k]))
        D=[[[C[i][j][k] for j in range(rows) if j!=k and j!=i] for k in range(rows)] for i in range(rows)]
        for i in range(rows):
            for k in range(rows):
                model.addConstr(sum(D[i][k])+B[i][k]>=1-M*(2-var_x[i]-var_x[j]))
                model.addConstr(sum(D[i][k])+B[i][k]<=1+M*(2-var_x[i]-var_x[j]))
        model.setObjective(sum(var_x), grb.GRB.MINIMIZE)
        model.optimize()
        VAR = [var_x[i].x for i in range(rows)]
        print(VAR)
        return VAR

    def CMDS_MARTIN_CP(self):
        rows = self.nodes
        M=self.M
        A = np.array(self.adj)
        A0 = np.array(self.adj) + np.array([[0 if j != i else 1 for j in range(rows)] for i in range(rows)])
        params = self.params
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 5 * 60)

        var_x=[]
        var_y=[]
        var_z=[]
        for i in range(rows):
            x=model.addVar(vtype=cp.COPT.BINARY, name='x'+str(i))
            var_x.append(x)
            var_1=[]
            var_2=[]
            for j in range(rows):
                y=model.addVar(vtype=cp.COPT.BINARY, name='y'+str(i)+str(j))
                var_1.append(y)
                var_3=[]
                for k in range(rows):
                    z=model.addVar(vtype=cp.COPT.BINARY, name='z'+str(i)+str(j)+str(k))
                    var_3.append(z)
                var_2.append(var_3)
            var_y.append(var_1)
            var_z.append(var_2)
        A0 = [[A[i][j] * var_x[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(sum(A0[i])>=1)
        B=[[A[i][j]*var_y[i][j] for j in range(rows)] for i in range(rows)]
        C=[[[A[i][j]*var_z[i][j][k] for k in range(rows)] for j in range(rows)] for i in range(rows)]
        model.addConstr(sum([sum(B[i]) for i in range(rows)])==sum(var_x)-1)
        for i in range(rows):
            for j in range(rows):
                model.addConstr(B[i][j]<=var_x[i])
                model.addConstr(B[i][j]<=var_x[j])
                for k in range(rows):
                    model.addConstr(C[i][j][k]<=B[i][j])
                    model.addConstr(C[i][j][k]<=var_x[k])
                    model.addConstr(C[j][i][k]<=B[i][j])
                    model.addConstr(C[j][i][k] <= var_x[k])
                    model.addConstr(B[i][j]-M*(3-var_x[i]-var_x[j]-var_x[k])-C[i][j][k]-C[j][i][k]<=0)
                    model.addConstr(C[i][j][k]+C[j][i][k]<=B[i][j]+M*(3-var_x[i]-var_x[j]-var_x[k]))
        D=[[[C[i][j][k] for j in range(rows) if j!=k and j!=i] for k in range(rows)] for i in range(rows)]
        for i in range(rows):
            for k in range(rows):
                model.addConstr(sum(D[i][k])+B[i][k]>=1-M*(2-var_x[i]-var_x[j]))
                model.addConstr(sum(D[i][k])+B[i][k]<=1+M*(2-var_x[i]-var_x[j]))
        model.setObjective(sum(var_x), cp.COPT.MINIMIZE)
        model.solve()
        VAR = [var_x[i].X for i in range(rows)]
        return VAR


            
    
    def CMDS_OR(self, cons):
        if cons=='MTZ':
            var=self.CMDS_MTZ_OR()
        elif cons=='SCF':
            var=self.CMDS_SCF_OR()
        return var
    
    def CMDS_GRB(self, cons):
        if cons=='MTZ':
            var=self.CMDS_MTZ_GRB()
        elif cons=='SCF':
            var=self.CMDS_SCF_GRB()
        elif cons=='MAR':
            var=self.CMDS_MARTIN_GRB()
        return var
    
    def CMDS_CP(self, cons):
        if cons=='MTZ':
            var=self.CMDS_MTZ_CP()
        elif cons=='SCF':
            print('Sorry! COPT can not solve this problem')
            return []
        elif cons=='MAR':
            var=self.CMDS_MARTIN_CP()
        return var
    
    #Solver MDS using the solver selectedby the user.
    def CMDS(self, cons):
        solver=self.solver
        if solver=='OR':
            var=self.CMDS_OR(cons)
        elif solver=='GRB':
            var=self.CMDS_GRB(cons)
        elif solver=='CP':
            var=self.CMDS_CP(cons)
        return var
    
    def PMDS_GRB(self,Z):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        
        var_1=[]
        for i in range(rows):
            x=model.addVar(vtype=grb.GRB.BINARY, name='x'+str(i))
            var_1.append(x)
        var_p=[]
        for i in range(rows):
            var_2=[]
            for j in range(rows):
                if j!=i:
                    p=model.addVar(vtype=grb.GRB.BINARY, name='p'+str(i)+str(j))
                    model.addConstr(p<=A[i][j])
                    var_2.append(p)
                else:
                    p=model.addVar(vtype=grb.GRB.BINARY, name='p'+str(i)+str(j))
                    model.addConstr(p==0)
                    var_2.append(p)
            var_p.append(var_2)
        B=A*var_1
        C=[[A[i][j]*var_p[i][j] for j in range(rows)] for i in range(rows)]
        Z0=A*Z
        Z1=Z0*np.array(var_p).T
        for i in range(rows):
            model.addConstr(sum(B[i])+sum(Z1[i])>=1)
            model.addConstr(sum(C[i])==Z[i])
        model.setObjective(sum(var_1), grb.GRB.MINIMIZE)
        model.optimize()
        VAR = [var_1[i].x for i in range(rows)]
        print(VAR)
        return VAR
    
    def PMDS_CP(self,Z):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        
        var_1=[]
        for i in range(rows):
            x=model.addVar(vtype=cp.COPT.BINARY, name='x'+str(i))
            var_1.append(x)
        var_p=[]
        for i in range(rows):
            var_2=[]
            for j in range(rows):
                if j!=i:
                    p=model.addVar(vtype=cp.COPT.BINARY, name='p'+str(i)+str(j))
                    model.addConstr(p<=A[i][j])
                    var_2.append(p)
                else:
                    p=model.addVar(vtype=cp.COPT.BINARY, name='p'+str(i)+str(j))
                    model.addConstr(p==0)
                    var_2.append(p)
            var_p.append(var_2)
        B=[[A[i][j]*var_1[j] for j in range(rows)] for i in range(rows)]
        C=[[A[i][j]*var_p[i][j] for j in range(rows)] for i in range(rows)]
        Z0=[[A[i][j]*Z[j] for j in range(rows)] for i in range(rows)]
        Z1=[[Z0[i][j]*var_p[j][i] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(sum(B[i])+sum(Z1[i])>=1)
            model.addConstr(sum(C[i])==Z[i])
        model.setObjective(sum(var_1), cp.COPT.MINIMIZE)
        status=model.solve()
        VAR = [var_1[i].x for i in range(rows)]
        if status==cp.COPT.OPTIMAL:
            print(VAR)
        return VAR
    
    def DRDP_1_GRB(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        
        var_x=[]
        var_y=[]
        var_z=[]
        for i in range(rows):
            x=model.addVar(vtype=grb.GRB.BINARY, name='x'+str(i))
            var_x.append(x)
            y=model.addVar(vtype=grb.GRB.BINARY, name='y'+str(i))
            var_y.append(y)
            z=model.addVar(vtype=grb.GRB.BINARY, name='z'+str(i))
            var_z.append(z)
        B=A*var_y
        C=A*var_z
        for i in range(rows):
            model.addConstr(var_x[i]+var_y[i]+var_z[i]+0.5*sum(B[i])+sum(C[i])>=1)
            model.addConstr(var_x[i]<=sum(B[i])+sum(C[i]))
            model.addConstr(var_x[i]+var_y[i]+var_z[i]<=1)
        model.setObjective(sum(var_x)+2*sum(var_y)+3*sum(var_z), grb.GRB.MINIMIZE)
        model.optimize()
        ans=[]
        for i in range(rows):
            if var_x[i].x==1:
                ans.append(1)
            elif var_y[i].x==1:
                ans.append(2)
            elif var_z[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_1_CP(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        
        var_x=[]
        var_y=[]
        var_z=[]
        for i in range(rows):
            x=model.addVar(vtype=cp.COPT.BINARY, name='x'+str(i))
            var_x.append(x)
            y=model.addVar(vtype=cp.COPT.BINARY, name='y'+str(i))
            var_y.append(y)
            z=model.addVar(vtype=cp.COPT.BINARY, name='z'+str(i))
            var_z.append(z)
        B=[[A[i][j]*var_y[j] for j in range(rows)] for i in range(rows)]
        C=[[A[i][j]*var_z[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(var_x[i]+var_y[i]+var_z[i]+0.5*sum(B[i])+sum(C[i])>=1)
            model.addConstr(var_x[i]<=sum(B[i])+sum(C[i]))
            model.addConstr(var_x[i]+var_y[i]+var_z[i]<=1)
        model.setObjective(sum(var_x)+2*sum(var_y)+3*sum(var_z), cp.COPT.MINIMIZE)
        model.solve()
        ans=[]
        for i in range(rows):
            if var_x[i].x==1:
                ans.append(1)
            elif var_y[i].x==1:
                ans.append(2)
            elif var_z[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_2_GRB(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        
        var_p=[]
        var_q=[]
        var_r=[]
        for i in range(rows):
            p=model.addVar(vtype=grb.GRB.BINARY, name='p'+str(i))
            var_p.append(p)
            q=model.addVar(vtype=grb.GRB.BINARY, name='q'+str(i))
            var_q.append(q)
            r=model.addVar(vtype=grb.GRB.BINARY, name='r'+str(i))
            var_r.append(r)
        B=A*var_q
        C=A*var_r
        for i in range(rows):
            model.addConstr(var_p[i]+0.5*sum(B[i])+0.5*sum(C[i])>=1)
            model.addConstr(var_q[i]+sum(B[i])>=var_p[i])
            model.addConstr(var_p[i]>=var_q[i])
            model.addConstr(var_q[i]>=var_r[i])
        model.setObjective(sum(var_p)+sum(var_q)+sum(var_r), grb.GRB.MINIMIZE)
        model.optimize()
        ans=[]
        for i in range(rows):
            if var_p[i].x==1 and var_q[i].x==0:
                ans.append(1)
            elif var_q[i].x==1 and var_p[i].x==1:
                ans.append(2)
            elif var_r[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_2_CP(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        
        var_p=[]
        var_q=[]
        var_r=[]
        for i in range(rows):
            p=model.addVar(vtype=cp.COPT.BINARY, name='p'+str(i))
            var_p.append(p)
            q=model.addVar(vtype=cp.COPT.BINARY, name='q'+str(i))
            var_q.append(q)
            r=model.addVar(vtype=cp.COPT.BINARY, name='r'+str(i))
            var_r.append(r)
        B=[[A[i][j]*var_q[j] for j in range(rows)] for i in range(rows)]
        C=[[A[i][j]*var_r[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(var_p[i]+0.5*sum(B[i])+0.5*sum(C[i])>=1)
            model.addConstr(var_q[i]+sum(B[i])>=var_p[i])
            model.addConstr(var_p[i]>=var_q[i])
            model.addConstr(var_q[i]>=var_r[i])
        model.setObjective(sum(var_p)+sum(var_q)+sum(var_r), cp.COPT.MINIMIZE)
        model.solve()
        ans=[]
        for i in range(rows):
            if var_p[i].x==1 and var_q[i].x==0:
                ans.append(1)
            elif var_q[i].x==1 and var_p[i].x==1:
                ans.append(2)
            elif var_r[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_3_GRB(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)
        
        var_y=[]
        var_z=[]
        for i in range(rows):
            y=model.addVar(vtype=grb.GRB.BINARY, name='y'+str(i))
            var_y.append(y)
            z=model.addVar(vtype=grb.GRB.BINARY, name='z'+str(i))
            var_z.append(z)
        B=A*var_y
        C=A*var_z
        for i in range(rows):
            model.addConstr(var_y[i]+var_z[i]+0.5*sum(B[i])+sum(C[i])>=1)
            model.addConstr(var_y[i]+var_z[i]<=1)
        model.setObjective(2*sum(var_y)+3*sum(var_z), grb.GRB.MINIMIZE)
        model.optimize()
        ans=[]
        for i in range(rows):
            if var_y[i].x==1:
                ans.append(2)
            elif var_z[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_3_CP(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        
        var_y=[]
        var_z=[]
        for i in range(rows):
            y=model.addVar(vtype=cp.COPT.BINARY, name='y'+str(i))
            var_y.append(y)
            z=model.addVar(vtype=cp.COPT.BINARY, name='z'+str(i))
            var_z.append(z)
        B=[[A[i][j]*var_y[j] for j in range(rows)] for i in range(rows)]
        C=[[A[i][j]*var_z[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(var_y[i]+var_z[i]+0.5*sum(B[i])+sum(C[i])>=1)
            model.addConstr(var_y[i]+var_z[i]<=1)
        model.setObjective(2*sum(var_y)+3*sum(var_z), cp.COPT.MINIMIZE)
        model.solve()
        ans=[]
        for i in range(rows):
            if var_y[i].x==1:
                ans.append(2)
            elif var_z[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_4_GRB(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        params=self.params
        env = grb.Env(params=params)

        # Create the model within the Gurobi environment
        model = grb.Model(env=env)
        model.setParam('TimeLimit', 30*60)

        var_q=[]
        var_r=[]
        for i in range(rows):
            q=model.addVar(vtype=grb.GRB.BINARY, name='q'+str(i))
            var_q.append(q)
            r=model.addVar(vtype=grb.GRB.BINARY, name='r'+str(i))
            var_r.append(r)
        B=A*var_q
        C=A*var_r
        for i in range(rows):
            model.addConstr(var_q[i]+0.5*sum(B[i])+0.5*sum(C[i])>=1)
            model.addConstr(var_q[i]>=var_r[i])
        model.setObjective(2*sum(var_q)+sum(var_r), grb.GRB.MINIMIZE)
        model.optimize()
        ans=[]
        for i in range(rows):
            if var_q[i].x==1 and var_r[i].x==0:
                ans.append(2)
            elif var_r[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
    
    def DRDP_4_CP(self):
        rows = self.nodes
        A = np.array(self.adj)
        A0=np.array(self.adj)+np.array([[0 if j!=i else 1 for j in range(rows)] for i in range(rows)])
        
        env = cp.Envr()

        # Create the model within the Gurobi environment
        model = env.createModel()
        model.setParam(cp.COPT.Param.TimeLimit, 30*60)
        
        var_q=[]
        var_r=[]
        for i in range(rows):
            q=model.addVar(vtype=cp.COPT.BINARY, name='q'+str(i))
            var_q.append(q)
            r=model.addVar(vtype=cp.COPT.BINARY, name='r'+str(i))
            var_r.append(r)
        B=[[A[i][j]*var_q[j] for j in range(rows)] for i in range(rows)]
        C=[[A[i][j]*var_r[j] for j in range(rows)] for i in range(rows)]
        for i in range(rows):
            model.addConstr(var_q[i]+0.5*sum(B[i])+0.5*sum(C[i])>=1)
            model.addConstr(var_q[i]>=var_r[i])
        model.setObjective(2*sum(var_q)+sum(var_r), cp.COPT.MINIMIZE)
        model.solve()
        ans=[]
        for i in range(rows):
            if var_q[i].x==1 and var_r[i].x==0:
                ans.append(2)
            elif var_r[i].x==1:
                ans.append(3)
            else:
                ans.append(0)
        return ans
        
    
    def PMDS(self, Z):
        if self.solver=='GRB':
            return self.PMDS_GRB(Z)
        elif self.solver=='CP':
            return self.PMDS_CP(Z)
    
    def DRDP(self, n):
        if self.solver=='GRB':
            if n==1:
                return self.DRDP_1_GRB()
            elif n==2:
                return self.DRDP_2_GRB()
            elif n==3:
                return self.DRDP_3_GRB()
            elif n==4:
                return self.DRDP_4_GRB()
        elif self.solver=='CP':
            if n==1:
                return self.DRDP_1_CP()
            elif n==2:
                return self.DRDP_2_CP()
            elif n==3:
                return self.DRDP_3_CP()
            elif n==4:
                return self.DRDP_4_CP()
                  
    
    #Visualize the result.
    def draw(self, var):
        A = np.array(self.adj)
        rows=self.nodes
        #Construct the graph using the adjacency matrix
        G = nx.from_numpy_array(A)

        # Optionally, you can set labels for the nodes
        # For example, if you want to label them 0, 1, 2, 3
        node_labels = {i: str(i) for i in range(A.shape[0])}
        
        dom=[i for i in range(rows) if var[i]==1]
        
        fig=plt.figure(figsize=(16,8))
        colors=['r' if i in dom else 'b' for i in range(rows)]
        nx.draw(G,labels=node_labels,with_labels=True, node_color=colors)
        plt.show()
        
    def draw_DRDP(self, var):
        A = np.array(self.adj)
        rows=self.nodes
        #Construct the graph using the adjacency matrix
        G = nx.from_numpy_array(A)

        # Optionally, you can set labels for the nodes
        # For example, if you want to label them 0, 1, 2, 3
        node_labels = {i: str(i) for i in range(A.shape[0])}
        colors=['r', 'b', 'y', 'g']
        col=[colors[var[i]] for i in range(rows)]
        
        
        fig=plt.figure(figsize=(16,8))
        nx.draw(G,labels=node_labels,with_labels=True, node_color=col)
        plt.show()