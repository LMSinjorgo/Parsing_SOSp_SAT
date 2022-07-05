# Parsing_SOSp_SAT
A MATLAB function for parsing the $SOS_p$ basis corresponding to a (MAX-)3-SAT instance

Given a logical 3-SAT proposition on $n$ variables $m$ clauses, the function $\texttt{SOS_p_Parse.m}$ takes as input matrix $Q \in \{ 0, 1 \}^{m \times 2n}$. This matrix satisfies, for $i \in [m]$ and $j \in [n]$:
$$Q_{i,j} = 1 \text{ if } x_j \text{ appears in clause } i$$
and for $i \in [m]$ and $n+1 \leq j \leq 2n$,
$$Q_{i,j} = 1 \text{ if } \neg x_{j-n} \text{ appears in clause } i.$$

The function $\texttt{cnfConverter.m}$ takes as input a $\texttt{.cnf}$ file, as specified by the DIMACS CNF file format (see https://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html), and converts this into a matrix $Q$ as described above.

The function $\texttt{createRandomInstance.m}$, with inputs $n$, $m$ and $k$ creates a $Q$ matrix corresponding to a (MAX-)SAT instance on $n$ variables and $m$ clauses, all of which have length $k$.
