# Integer-Programming-Dominating-Set
Implement integer programming for the minimum dominating set problem and its variants using CP_Solver, Gurobi, and COPT.

1.Introduction
In graph theory, a dominating set for a graph $G$ is a subset $D$ of its vertices, such that any vertex of $G$ is in $D$, or has a neighbor in $D$. The minimum dominating set problem concerns the size of the smallest dominating set of the input graph $G$, denoted by $\gamma(G)$. The application of domination in graph lies in various fields in solving real life problems. It includes social network theory, land surveying, radio stations, computer communication networks, school bus routing, sets of representatives, interconnection networks, etc. Depends on the specific application, domination problem has various kinds of variants, like connected domination, strong domination, and Roman domination, etc.

Integer programming is an important tool for solving the domination problem. With the help of optimization solvers, we can solve the dominating set problem efficiently. In this project, we model some variants of dominating set using integer programming and solve them using CP SAT Solver from OR\_TOOLS, Gurobi and COPT.
