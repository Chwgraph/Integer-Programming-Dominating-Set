# Integer-Programming-Dominating-Set
Implement integer programming for the minimum dominating set problem and its variants using CP_Solver, Gurobi, and COPT.

# 1.Introduction

In graph theory, a dominating set for a graph $G$ is a subset $D$ of its vertices, such that any vertex of $G$ is in $D$, or has a neighbor in $D$. The minimum dominating set problem concerns the size of the smallest dominating set of the input graph $G$, denoted by $\gamma(G)$. The application of domination in graph lies in various fields in solving real life problems. It includes social network theory, land surveying, radio stations, computer communication networks, school bus routing, sets of representatives, interconnection networks, etc. Depends on the specific application, domination problem has various kinds of variants, like connected domination, strong domination, and Roman domination, etc.

Integer programming is an important tool for solving the domination problem. With the help of optimization solvers, we can solve the dominating set problem efficiently. In this project, we model some variants of dominating set using integer programming and solve them using CP SAT Solver from OR\_TOOLS, Gurobi and COPT.

# 2.Licenses Required

For Gurobi, a WLS license is needed. Before compiling the program, please make sure you input your WLSACCESSID, WLSSECRECT and LICENSEID in the params. For COPT, a free trial license for academic or business use is needed. See the application link here: [Application Link](https://www.shanshu.ai/copt). CP SAT Solver is free, but needed to install the library 'ortools' first.

# 3. Variants of Domination Problem Covered
   
This project covers three variants of domination problem in total. We briefly introduce how each domination problem is formulated in our program here. For details about each constraint, see the references [5] and [6].

## A) Minimum Dominating Set(MDS) Problem

For a graph $G$, a vertex subset $D$ is called a dominating set of $G$ if for any vertex $v$ in $V(G)$, either it has a neighbour in $D$ or it is in $D$ itself. The minimum dominating set problem is to find a dominating set with smallest cardinality. The decision version of the MDS is a classical NP-complete problem. [1] An integer programming formation for MDS problem is as follows: 

min $\sum_{v \in V(G)} x_v$

s.t. $\sum_{u \in N(v)} x_u + x_v >=1, \forall v \in V(G)$

MDS problem can be easily generalized by k-MDS problem in which we require each vertex to have at least $k$ neighbours in the dominating set. This is called k-MDS Prolem.

## B) Minimum Strongly/Total Dominating Set (MSDS) Problem

A total dominating set (or strongly-dominating set) is a set of vertices such that all vertices in the graph, including the vertices in the dominating set themselves, have a neighbor in the dominating set.[2] An integer programming formation for SMDS problem is as follows:

min $\sum_{v \in V(G)} x_v$

s.t. $\sum_{u \in N(v)} x_u >=1, \forall v \in V(G)$

## C) Minimum Connected Dominating Set (MCDS) Problem

MCDS problem requires the minimum dominating set obtained can induce a connected subgraph. There are various kinds of connectivity constraints. 

### 1) Miller-Tucker-Zemlin (MTZ) Constraints

The MTZ constraint was first proposed for the travel salesman problem [3], and were used to eliminate sub-tours when solving the kcardinality tree problem [4]. It is formulated as below:

$\sum_{i \in V} y_{n+2, i}=1$ 

$\sum_{i:(i,j) \in A} y_{ij}=1, \forall v \in V, A= \{ (n+1, n+2)\} \cup \{\cup_{i=1}^{n} \{(n+1, i), (n+2, i)\}\} \cup E \cup E', E'=\{(j,i)| (i,j) \in E\, i>j\}$

$y_{n+1, i}+y_{i,j} \leq 1, \forall (i,j) \in E \cup E'$

$(n+1)y_{ij}+u_i-u_j+(n-1)y_{ji} \leq n, \forall (i,j) \in E \cup E'$

$(n+1)y_{ij}+u_i-u_j \leq n, \forall (i,j) \in A-(E \cup E')$

$y_{n+1, n+2}=1$

$u_{n+1}=0$

$1 \leq u_i \leq n+1, \forall i \in V \cup \{n+2\}$

$x_i=1-y_{n+1, i}, \forall i \in V$

### 2) Martin Constraints

In [7], Martin presented a reformulation for solving the minimum spanning tree problem with a polynomial number of constraints instead of an exponential number of constraints.

$\sum_{(i,j) \in E}=\sum_{i \in V} x_i -1$

$y_{ij} \leq x_i, y_{ij} \leq x_j, \forall (i,j) \in E$

$z_{ij}^{k} \leq y_{ij}, z_{ij}^{k} \leq x_k, \forall (i,j) \in E, k \in V$

$z_{ji}^{k} \leq y_{ij}, z_{ji}^{k} \leq x_k, \forall (i,j) \in E, k \in V$

$y_{ij}-M(3-x_i-x_j-x_k) \leq z_{ij}^{k}+z_{ji}^{k} \leq y_{ij} +M(3-x_i-x_j-x_k), \forall i, j, k \in V$

$1-M(2-x_i-x_j) \leq \sum_{k' \in V-\{i,j\}} z_{ik'}^{j} +y_{ij} \leq 1+M(2-x_i-x_j), \forall i,j \in V$

$y_{ij}, z_{ij}^{k} \in \{0,1\}, \forall (i,j) \in E$

$y_{ij}=0, z_{ij}^{k}=0, \forall i,j,k \in V, (i,j) \notin E$

### 3) Single-Commodity Flow (SCF) Constraints

We define E' as we did in 1).

$\sum_{i \in V} r_i =1$

$r_i \leq x_i, \forall i \in V$

$f_{ij} \geq 0, \forall (i,j) \in E \cup E'$

$f_{ij} \leq x_i \sum_{k \in V} x_k, f_{ij} \leq x_j \sum_{k \in V} x_k, \forall (i,j) \in E \cup E'$

$\sum_{j} f_{ij} \leq n(1-r_i), \forall v \in V$

$\sum_{j} f_{ji} -\sum_{i} f_{ij}=x_i-r_i \sum_{j} x_j, \forall i \in V$

$r_i \in \{0,1\}, \forall v \in V$

## D) Minimum Power Dominating Set (MPDS) Problem

For a power graph $G = (V,E)$, there is a given subset $V_Z \subset V$ of zero-injection vertices. As explained in [8], a power dominating set $D$ is obtained by leveraging two physical laws: (1) if $v \in D$, then v and its neighbors (denoted by $N(v)$) are all covered (Ohm’s law); (2) if $v ∈ V_Z$, and all vertices within the set $\{v\} \cup N(v)$ except one are covered, then the uncovered vertex in $\{v\} \cup N(v)$ is also covered (Kirchhoff’s current law). The power dominating set (PDS) problem is to find a subset $D$ with smallest cardinality that covers all vertices in $V$. In [9], a formulation for MPDS problem is proposed as follows:

$\min_{x,p} \sum_{i} x_i$

$s.t. \sum_{j} a_{ij}x_j + \sum_{j} a_{ij}Z_jp_{ij} \geq 1, \forall i \in V$

$\sum_{j} a_{ij}p_{ij} =Z_i, \forall i \in V$

$p_{ij}=0, \forall i,j \quad \text{if} \quad a_{ij}=0 \quad \text{or} \quad i \notin V_Z$

$x_i, p_{ij} \in \{0,1\}, \forall i, j \in V$

## E) Double Roman Domination Problem (DRDP)

A double Roman dominating function (DRDF) is a function $f : V \longrightarrow \{0, 1, 2, 3\}$ having the property that if $f(v) = 0$, then vertex $v$ must have at least two neighbors assigned $2$ under $f$ or at least one neighbor $u$ with $f(u) = 3$, and if $f(v) = 1$, then vertex $v$ must have at least one neighbor $u$ with $f(u) \geq 2$. The weight of a DRDF $f$ is the sum $w(f) = \sum_{v} f(v)$, and the double Roman domination number of a graph $G$, denoted by $\gamma_{dR}(G)$, is the minimum weight among all possible DRDFs.

In this project, we use four formulations proposed in [6] for DRDP.

### 1) DRDP-1

For each vertex $v \in V$, we define $x_v$ as the indicator variable for $f(v)=1$, $y_v$ as the indicator variable for $f(v)=2$, $z_v$ as the indicator variable for $f(v)=3$.

$\min \sum_{v} x_v + \sum_{v} 2y_v + \sum 3 z_v$

$s.t. x_v+y_v+z_v+\frac{1}{2}\sum_{u \in N(v)} y_u + \sum_{u \in N(v)} z_u \geq 1, \forall v \in V$

$\sum_{u \in N(v)} (y_u+z_u) \geq x_v, \forall v \in V$

$x_v+y_v+z_v \leq 1, \forall v \in V$

$x_v, y_v, z_v \in \{0,1\}, \forall v \in V$


### 2) DRDP-2

For each vertex $v \in V$, we define $p_v$ as the indicator variable for $f(v) \geq 1$, $q_v$ as the indicator variable for $f(v) \geq 2$, $r_v$ as the indicator variable for $f(v)=3$.

$\min \sum_{v} p_v + \sum_{v} q_v + \sum r_v$

$s.t. p_v+\frac{1}{2}\sum_{u \in N(v)} q_u + \sum_{u \in N(v)} r_u \geq 1, \forall v \in V$

$q_v+ \sum_{u \in N(v)} q_u \geq p_v, \forall v \in V$

$r_v \leq q_v \leq p_v, \forall v \in V$

$r_v, p_v, q_v \in \{0,1\}, \forall v \in V$


### 3) DRDP-1'

Improved version of DRDP-1.

$\min \sum_{v} 2y_v + \sum 3 z_v$

$s.t. y_v+z_v+\frac{1}{2}\sum_{u \in N(v)} y_u + \sum_{u \in N(v)} z_u \geq 1, \forall v \in V$

$\sum_{u \in N(v)} (y_u+z_u) \geq x_v, \forall v \in V$

$y_v+z_v \leq 1, \forall v \in V$

$y_v, z_v \in \{0,1\}, \forall v \in V$


### 4) DRDP-2'

Improved version of DRDP-2.

$\min \sum_{v} 2q_v + \sum r_v$

$s.t. q_v+\frac{1}{2}\sum_{u \in N(v)} q_u + \frac{1}{2} \sum_{u \in N(v)} r_u \geq 1, \forall v \in V$

$r_v \leq q_v, \forall v \in V$

$r_v, q_v \in \{0,1\}, \forall v \in V$


# 4. Guidance to Use Program

## A) Basic Parameters

|Name|Description|
|---|---|
|adj_matrix|The adjacency matrix of your input graph|
|nodes|Number of vertices in the graph|
|solver|The solver you would like to use: 'OR', 'GRB', or 'CP'|
|params|The params for license of Gurobi if you choose 'GRB'|
|M|Value of big M, defaulted value is 10000|

## B) Functions

|Name|Description|Parameters|
|---|---|---|
|MDS|Function for solving strong/basic k-MDS problem|k: value for k in k-MDS problem, ver: 's' for strong, 'w' for basic|
|CMDS|Function for solving CMDS Problem|cons: The connected constraints you would like to use, 'MTZ' for Miller-Tucker-Zemlin, 'MAR' for Martin, and 'SCF' for Single-Commodity Flow|
|PMDS|Function for solving PMDS Problem|Z:zero-injection vertex set indicated using a 0-1 Numpy Array|
|DRDP|Function for solving DRDP|n: No. of the constraints you would like to use, 1 for DRDP-1, 2 for DRDP-2, 3 for DRDP-1' and 4 for DRDP-2'|
|draw|Function to visualize the result by drawing the graph with vertices in the dominating set colored red|var: vector output by one of the functions: MDS, CMDS and PMDS|
|draw_DRDP|Function to visualize the result of DRDP by drawing the graph with vertices f(v)=0 red, f(v)=1 blue, f(v)=2 yellow and f(v)=3 green|var: vector output by DRDP|

## C) Notes

1. For PDMS and Martin Constraints for CDMS, only Gurobi and COPT are available;
2. For SCF Constraints for CDMS, only Gurobi and CP SAT Solver are available;
3. The time limit for each problem is set 5 minutes.
4. For Gurobi params, please input your params as suggested on the official website of Gurobi;
5. Comments and suggestions are welcomed. Since this is the original version, there may be some issues.

# References
1. Garey, M.R., Johnson, D.S.: Computers and Intractability: A Guide to the Theory of NPCompleteness. W. H. Freeman & Co., New York (1979)
2. West, Douglas B. (2001), Introduction to Graph Theory (2 ed.), Pearson Education.
3. Miller, C.E., Tucker, A.W., Zemlin, R.A.: Integer programming formulation of traveling salesman problems. J. Assoc. Comp. Mach. 7, 326–329 (1960)
4. Quintao, F.R., da Cunha, A.S., Mateus, G.R., Lucena, A.: The k-cardinality tree problem: reformulations and lagrangian relaxation. Disc. Appl. Math. 158, 1305–1314 (2010)
5. Fan, Neng, and Jean-Paul Watson. "Solving the connected dominating set problem and power dominating set problem by integer programming." Combinatorial Optimization and Applications: 6th International Conference, COCOA 2012, Banff, AB, Canada, August 5-9, 2012. Proceedings 6. Springer Berlin Heidelberg, 2012.
6. Cai, Qingqiong, et al. "Integer linear programming formulations for double roman domination problem." Optimization Methods and Software 37.1 (2022): 1-22.
7. Martin, R.K.: Using separation algorithms to generate mixed integer model reformulations. Oper. Res. Lett. 10, 119–128 (1991)
8. Aazami, A.: Domination in graphs with bounded progagation: algorithms, formulations and hardness results. J. Comb. Optim. 19, 429–456 (2010)
9. Aminifar, F., Khodaei, A., Fotuhi-Firuzabad, M., Shahidehpour, M.: Contingencyconstrained PMU placement in power networks. IEEE Trans. Power Syst. 25, 516–523 (2010)
