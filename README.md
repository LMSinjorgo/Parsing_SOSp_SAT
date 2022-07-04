# Parsing_SOSp_SAT
A MATLAB function for parsing the $SOS_p$ basis corresponding to a (MAX-)3-SAT instance

Given a logical 3-SAT proposition on $n$ variables $m$ clauses, this function takes as input matrix $Q \in \{ 0, 1 \}^{m \times 2n}$. This matrix satisfies, for $i \in [m]$ and $j \in [n]$:
$$Q_{i,j} = 1 \text{ if } x_j \text{ appears in clause } i$$
and for $i \in [m]$ and $n+1 \leq j \leq 2n$,
$$Q_{i,j} = 1 \text{ if } \neg x_{j-n} \text{ appears in clause } i.$$



