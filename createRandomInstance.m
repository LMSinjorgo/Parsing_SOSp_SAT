function [Q] = createRandomInstance(n,m, k)
% creates a random SAT problem instance
% with n variables and m clauses

% the SAT instances will only have clauses, of which the length
% is contained in clauseLengths

clauses = [];
for j = 1:m
    clauses = [clauses; randsample(n,k)];
end
bernouilliDraw = binornd(1,0.5,m*k,1); 
clauses(bernouilliDraw == 1) = clauses(bernouilliDraw == 1)+n;
i = repelem(1:m,k)';
v = ones(m*k,1);
Q = sparse(i, clauses,v ,m,2*n);


end

