# Integer-Programming-Dominating-Set
Implement integer programming for the minimum dominating set problem and its variants using CP_Solver, Gurobi, and COPT.

1.Introduction

In graph theory, a dominating set for a graph $G$ is a subset $D$ of its vertices, such that any vertex of $G$ is in $D$, or has a neighbor in $D$. The minimum dominating set problem concerns the size of the smallest dominating set of the input graph $G$, denoted by $\gamma(G)$. The application of domination in graph lies in various fields in solving real life problems. It includes social network theory, land surveying, radio stations, computer communication networks, school bus routing, sets of representatives, interconnection networks, etc. Depends on the specific application, domination problem has various kinds of variants, like connected domination, strong domination, and Roman domination, etc.

Integer programming is an important tool for solving the domination problem. With the help of optimization solvers, we can solve the dominating set problem efficiently. In this project, we model some variants of dominating set using integer programming and solve them using CP SAT Solver from OR\_TOOLS, Gurobi and COPT.

2.Licenses Required

For Gurobi, a WLS license is needed. Before compiling the program, please make sure you input your WLSACCESSID, WLSSECRECT and LICENSEID in the params. For COPT, a free trial license for academic or business use is needed. See the application link here: https://www.shanshu.ai/copt. CP SAT Solver is free, but needed to install the library 'ortools' first.

3. Variants of Domination Problem Covered
   
This project covers three variants of domination problem in total.

A) Minimum Domination Set(MDS) Problem

For a graph $G$, a vertex subset $D$ is called a dominating set of $G$ if for any vertex $v$ in $V(G)$, either it has a neighbour in $D$ or it is in $D$ itself. The minimum dominating set problem is to find a dominating set with smallest cardinality. The decision version of the MDS is a classical NP-complete problem. [1] An integer programming formation for MDS problem is as follows: 

min $\sum_{v \in V(G)} x_v$

s.t. $\sum_{u \in N(v)} x_u + x_v >=1, \forall v \in V(G)$

References
1. Garey, M.R., Johnson, D.S.: Computers and Intractability: A Guide to the Theory of NPCompleteness. W. H. Freeman & Co., New York (1979)
