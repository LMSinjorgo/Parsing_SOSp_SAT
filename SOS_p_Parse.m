function ADMMinput = SOS_p_Parse(Q)

%%%%% INPUT

% Q matrix, see README file

%%%%% OUTPUT

%%% matrixSize: the length of monomial basis x.

%%% constant: the constant term of F_b

%%% unitMonomials: matrix with 4 columns. The nonzero entries in each row describe the
% corresponding monomial. these are the monomials that appear twice in the matrix
% X = x* x', where x is the monomial basis. The monomials in unitMonomials
% thus correspond to unit constraints.

% unitIdx: the linear indices corresponding to the unitMonomials. Each
% unitMonomial constitutes 1 index, i.e. unitIdx and unitMonomials
% have equal number of rows.

%%%% nonUnitMonomials: matrix with 4 columns. The nonzero entries in each row describe the
% corresponding monomial. The monomials in this matrix appear strictly more
% than 2 times in the matrix X = x*x', where x is the monomial basis.

%%%% paddedNonUnitIdx: each row contains the linear indices corresponding
% to the monomials in nonUnitMonomials. To make sure every row has equal
% length, the rows are padded with 1's.

%%%% coefficients: the coefficients of F_b corresponding the monomials
% in nonUnitMonomials.


n = size(Q,2)/2;
m = size(Q,1);
A = Q(:,1:n) - Q(:,(n+1):2*n);

%%%%%% compute F_B
ell = full(sum(abs(A),2));
constant = sum(0.5 .^ ell);
quadCounter = zeros(n,n);
linearCounter = zeros(n,1);
%%% check all the unique length 3 clauses
uniqueVariableClauses = unique(abs(A),'rows');
ellUnique = sum(uniqueVariableClauses,2);

unique3clauses = uniqueVariableClauses(ellUnique == 3,:);
numUniq3Clauses = size(unique3clauses,1);

% cubicCounter checks all the possible cubics appearing in F_B
% the first 3 columns refer to the indices of the monomials
% the last column takes the coefficient
cubicCounter = find(unique3clauses') - repelem(0:n:n*numUniq3Clauses-1,3)';
cubicCounter = reshape(cubicCounter,3,numUniq3Clauses)';
cubicCounter = [cubicCounter, zeros(numUniq3Clauses,1)];
for j = 1:m
    aj = A(j,:);
    vars = find(aj ~= 0);
    linearCounter(vars) = linearCounter(vars) - aj(vars)'/(2^ell(j));
    
    if ell(j) == 2
        quadCounter(vars(1),vars(2)) = quadCounter(vars(1),vars(2)) + aj(vars(1))*aj(vars(2)) / 4;
    elseif ell(j) == 3
        % divide by 8 for (1/2^3) = 1/8
        quadCounter(vars(1),vars(2)) = quadCounter(vars(1),vars(2)) + aj(vars(1))*aj(vars(2)) / 8;
        quadCounter(vars(1),vars(3)) = quadCounter(vars(1),vars(3)) + aj(vars(1))*aj(vars(3)) / 8;
        quadCounter(vars(2),vars(3)) = quadCounter(vars(2),vars(3)) + aj(vars(2))*aj(vars(3)) / 8;
        [~,matchingIdx] = ismember(vars, cubicCounter(:,1:3) ,'rows');
        cubicCounter(matchingIdx,4) = cubicCounter(matchingIdx,4) - aj(vars(1))*aj(vars(2))*aj(vars(3)) / 8;
    end
    
end



% create the monomial basis x
x = sparse(0,n);
matrixSize = 0;
for j = 1:size(uniqueVariableClauses,1)
    if ellUnique(j) == 2
        x = [x; uniqueVariableClauses(j,:)];
        matrixSize = 1+matrixSize;
    end
    if ellUnique(j) == 3
        vars = find(uniqueVariableClauses(j,:) == 1);
        quads = nchoosek(vars,2);
        for k = 1:3
            x(matrixSize+k,quads(k,:)) = 1;
        end
        matrixSize = matrixSize + 3;
    end
end

x = unique(x,'rows');
currentSize = size(x,1);

minSize = 1650 - n-1;
%minSize = -1;
deg = 2;
if currentSize < minSize
   % add all 2-degree combinations of most occuring vars
   varOccurence = sum(abs(A),1);
   k = ceil(sqrt( 2*(minSize - currentSize) ));
   % currentSize + binom(k,2) \approx minSize
   [~, maxOccuringVars] = maxk(varOccurence,k);
   combs = nchoosek(maxOccuringVars,deg);
   numCombs = size(combs,1);
   combs = combs'; combs = combs(:);
   
   
   % create sparse matrix
   rowIdx = repelem(1:numCombs,deg);
   newX = sparse(rowIdx',combs,ones(1,numCombs*2),numCombs,n);
   x = [x; newX];
end

x = unique(x,'rows');

x = [zeros(1,n); speye(n); x];
matrixSize = size(x,1);

% form the singleton idx,
% these indices keep track which elements of basis x contain a given
% singleton
singletonGroupSizes = full(sum(x,1));
idxHelp = repmat(1:matrixSize,1,n);
singletonIdx = idxHelp( x(:) == 1);

singletonMonomials = (1:n)';


%%%% end singletonIdx creation



% compute all the cross products between monomials
C1 = repelem(x(2:end,:),1:matrixSize-1,1);

numCrossPairs = size(C1,1);

s = x * ones(n,1);
numNonZeros =  ((size(x,1):-1:1)-1) * s;
s = cumsum(s); s(1) = 0; sIdx = cumsum(s);
[c,r] = find(x');
sparseRow = zeros(numNonZeros,1);
sparseCol = zeros(numNonZeros,1);

rowOffset = 1;
for j = 2:matrixSize-1
    sparseRow(sIdx(j-1)+1:sIdx(j)) = r(1:s(j))+rowOffset;
    sparseCol(sIdx(j-1)+1:sIdx(j)) = c(1:s(j));
    rowOffset = rowOffset+j;
end
C2 = sparse(sparseRow,sparseCol, ones(numNonZeros,1), numCrossPairs,n);

crossPairs = mod(C1+C2,2);

monomialSizes = full(sum(crossPairs,2));

upperTriuIdx = find(triu(ones(matrixSize),1))';

allIdx = cell(numCrossPairs,1);

groupSizes = zeros(numCrossPairs,1);
% groupSizes keeps track of the number of indices per monomial
if n > 4
    allMonomials = zeros(nchoosek(n,4)*4 + nchoosek(4,3)*3+nchoosek(4,2)*2,1);
else
    allMonomials = zeros(100,1);
end


varsCounter = 1;
idxCounter = 0;
degCount = zeros(4,1);

% for k = 1
monSizeLogic = monomialSizes == 1;
currentPairs = crossPairs(monSizeLogic,:);
currentIdx = upperTriuIdx(monSizeLogic);

for p = 1:n
    matchingIdx = currentPairs(:,p) == 1;
    allMonomials(varsCounter) = p;
    varsCounter = varsCounter+1;
    
    idxCounter = idxCounter + 1;

    allIdx{idxCounter,1} = currentIdx(matchingIdx);
    groupSizes(idxCounter) = size(allIdx{idxCounter,1},2);
end
degCount(1) = n;

monomialSizes(monSizeLogic) = [];
crossPairs(monSizeLogic,:) = [];
upperTriuIdx(monSizeLogic) = [];

for k = 2:3
    monSizeLogic = monomialSizes == k;
    currentPairsBig = crossPairs(monSizeLogic,:);
    currentIdxBig = upperTriuIdx(monSizeLogic);
    
    vectorOnes = ones(k-1,1);
    computeStep = k-1;
    for p = 1:n
        matchingVars = currentPairsBig(:,p) == 1;
        currentIdx  = currentIdxBig(matchingVars);
        
        currentPairs = currentPairsBig(matchingVars,:);
        numCurrentPairs = size(currentPairs,1);
        
        % find the variables per clause
        allVars = find(currentPairs') - repelem(0:n:n*(numCurrentPairs-1),k)';
        allVars = reshape(allVars,k,numCurrentPairs)';
        allVars = allVars(:,2:k);
        indexFind = (1:numCurrentPairs)';
        
        entryChecked = zeros(numCurrentPairs,1);
        while ~isempty(currentPairs)
            for j = 1:min(10000,size(currentPairs,1))
                if entryChecked(j) == 0
                    vars = allVars(j,:);
                    
                    products = currentPairs(:,vars)*vectorOnes;
                    
                    matchingIdx = indexFind(products == computeStep);
                    
                    allMonomials(varsCounter) = p;
                    allMonomials(varsCounter+1:varsCounter+computeStep) = vars;
                    
                    varsCounter = varsCounter + k;
                    
                    idxCounter = idxCounter + 1;

                    allIdx{idxCounter,1} = currentIdx(matchingIdx);
                    groupSizes(idxCounter) = size(matchingIdx,1);
                    % indicate the monomials that have just been
                    % placed in a group.
                    entryChecked(matchingIdx) = 1;
                    
                    %%% perhaps an idea to remove entries of currentPairs
                    % every set number of iterations???
                end
            end
            if j ~= size(currentPairs,1)
                idxHelp = 1:size(currentPairs,1);
                checkedEntries = idxHelp(entryChecked == 1);
                entryChecked(checkedEntries,:) = [];
                currentPairs(checkedEntries,:)=[];
                currentIdx(checkedEntries) = [];
                allVars(checkedEntries,:) = [];
            else
                currentPairs = [];
            end
        end
        currentPairsBig(matchingVars,:) = [];
        currentIdxBig(matchingVars) = [];
    end
    degCount(k) = idxCounter;
    monomialSizes(monSizeLogic) = [];
    crossPairs(monSizeLogic,:) = [];
    upperTriuIdx(monSizeLogic) = [];
end
% big matching loop end, k = 2:3

for k = 4
    % no need to first select only monomials with degree 4. All other
    % degrees have already been deleted from crossPairs and upperTriuIdx
    vectorOnes = ones(k-2,1);
    computeStep = k-1;
    for p = 1:(n-1)
        matchingVars = crossPairs(:,p) == 1;
        start = find(matchingVars == 1,1,'first');
        endPoint = start + nnz(matchingVars) - 1;
        currentIdx  = upperTriuIdx(start:endPoint);
        if isempty(currentIdx)
            continue;
        end
        currentPairs = crossPairs(start:endPoint,:);

        for q = p+1:n
            matchingVarsSmall = currentPairs(:,q) == 1;
            currentIdxSmall = currentIdx(matchingVarsSmall);
            if isempty(currentIdxSmall)
                continue;
            end
            currentPairsSmall = currentPairs(matchingVarsSmall,:);
            
            
            numCurrentPairs = size(currentPairsSmall,1);
            
            allVars = find(currentPairsSmall') - repelem(0:n:n*(numCurrentPairs-1),k)';
            allVars = reshape(allVars,k,numCurrentPairs)';
            allVars = allVars(:,3:k);
            
            entryChecked = zeros(numCurrentPairs,1);
            indexFind = (1:numCurrentPairs)';
            for j = 1:numCurrentPairs
                if entryChecked(j) == 0
                    vars = allVars(j,:);
                    
                    products = currentPairsSmall(:,vars)*vectorOnes;
                    
                    matchingIdx = indexFind(products == 2);

                    allMonomials(varsCounter) = p;
                    allMonomials(varsCounter+1) = q;
                    allMonomials(varsCounter+2:varsCounter+computeStep) = vars;
                    
                    varsCounter = varsCounter + k;
                    
                    idxCounter = idxCounter + 1;

                    allIdx{idxCounter,1} = currentIdxSmall(matchingIdx);
                    groupSizes(idxCounter) = size(matchingIdx,1);
                    
                    % indicate the monomials that have just been
                    % placed in a group.
                    entryChecked(matchingIdx) = 1;
                    
                    %%% perhaps an idea to remove entries of currentPairs
                    % every set number of iterations???
                end
            end
            currentPairs(matchingVarsSmall,:) = [];
            currentIdx(matchingVarsSmall) = [];
        end
        crossPairs(matchingVars,:) = [];
        upperTriuIdx(matchingVars) = [];
    end
    degCount(4) = idxCounter;
end


allMonomials = allMonomials(1:varsCounter-1);
idxVector = [1:degCount(1), repelem(degCount(1)+1:degCount(2),2),repelem(degCount(2)+1:degCount(3),3)];
if degCount(3)+1 <= degCount(4)
    idxVector = [idxVector, repelem(degCount(3)+1:degCount(4),4)];
end
allMonomials = sparse(idxVector, allMonomials',ones(size(allMonomials)),idxCounter,n);


monomialsDegrees = repelem((1:4)',[degCount(1), degCount(2)-degCount(1), degCount(3)-degCount(2),degCount(4)-degCount(3)]);
allIdx = allIdx(1:idxCounter,1);
groupSizes = groupSizes(1:idxCounter);

% only store the indices, not the coefficients since these are all 0
% anyways.


helpCell = allIdx(groupSizes == 1,1);
unitIdx= cat(1,helpCell{:});

allIdx = allIdx(groupSizes ~= 1,1);
nonPaddedIdx = horzcat(allIdx{:});


% get the nonunit Idx monomials
% non unit, as in, they appear strictly more than 2 times in the symmetric
% matrix x*x'. Their indices are given by paddedNonUnitIdx
nonUnitMonomials = allMonomials(groupSizes ~= 1,:);
nonUnitDegrees = monomialsDegrees(groupSizes ~= 1);

allMonomials(groupSizes ~= 1,:) = [];
unitMonomials = allMonomials;
unitDegrees = monomialsDegrees(groupSizes == 1);


% prepare for output
groupSizes = groupSizes(groupSizes ~= 1);


% we match the coefficients of F_b to allMonomials


% match the terms corresponding to quadratics, i.e.,
% x_1*x_2, x_3*x_5. Determine here already to determine the size of coeffs
[row,col] = find(quadCounter);
quadIdx = [row,col]; quadCounter = quadCounter(quadCounter ~= 0) / 2;
monomialSizes = nonUnitMonomials * ones(n,1);

coeffs = zeros(n+size(quadIdx,1)+size(cubicCounter,1),2);

% match the terms corresponding to singletons, i.e., x_1, x_2, x_3 etc
% divide by 2 for easier implementation
linearCounter = linearCounter/2;
for j = 1:n
    matchingIdx = find(nonUnitMonomials(:,j) == 1,1,'first');
    coeffs(j,:) = [matchingIdx, linearCounter(j)];
end


quadMonomials =  nonUnitMonomials(monomialSizes == 2,:);
cubicMonomials = nonUnitMonomials(monomialSizes == 3,:);

numQuadMonomials = size(quadMonomials,1);

for j = 1:size(quadIdx,1)
    b = quadMonomials(:,quadIdx(j,:) ) * ones(2,1);
    matchingIdx = find(b == 2,1,'first') + n;
    coeffs(n+j,:) = [matchingIdx, quadCounter(j)];
end

% now match the cubic terms to allMonomials
% divide by 2 for easier implementation
cubicCounter(:,4) = cubicCounter(:,4)/2;
numPrevMonomials = n + numQuadMonomials;
coeffIdxStart = n+size(quadIdx,1);
for j = 1:size(cubicCounter,1)
    b = cubicMonomials( :,cubicCounter(j,1:3)) * ones(3,1);
    matchingIdx = find(b == 3,1,'first') + numPrevMonomials;
    coeffs(coeffIdxStart+j,:) = [matchingIdx, cubicCounter(j,4)];
end

totalCoeffsSize = size(coeffs,1);

% coeffIndicator will be a boolean vector
% a 1 indicates the the corresponding monomial could possibly
% have a nonzero coefficient. This is useful for branching,
% where the coefficients have to be recomputed and reassigned.
coeffIndicator = sparse(coeffs(:,1),1,1,size(nonUnitMonomials,1),1,totalCoeffsSize + numQuadMonomials - size(quadIdx,1));
coeffIndicator = coeffIndicator + (monomialSizes == 2);
coeffIndicator = coeffIndicator ~= 0;

coeffs = sparse(coeffs(:,1),1,coeffs(:,2),size(nonUnitMonomials,1),1);


% prepare output for ADMM function
ADMMinput.unitMonomials = unitMonomials;
ADMMinput.unitIdx = unitIdx;
ADMMinput.unitDegrees = unitDegrees;

ADMMinput.nonUnitMonomials = nonUnitMonomials;
ADMMinput.nonUnitDegrees = nonUnitDegrees;
ADMMinput.groupSizes = groupSizes;
ADMMinput.nonPaddedIdx = nonPaddedIdx;


ADMMinput.coeffs = coeffs;
ADMMinput.coeffIndicator = coeffIndicator;
ADMMinput.x = x;
ADMMinput.numClauses = m;
ADMMinput.numVariables = n;
ADMMinput.constant = constant;


ADMMinput.singletonIdx = singletonIdx;
ADMMinput.singletonGroupSizes = singletonGroupSizes;
ADMMinput.singletonMonomials = singletonMonomials;
% small unit test
if size(unitIdx,1) + size(nonPaddedIdx,2) ~= matrixSize*(matrixSize-1) / 2
    error("not all indices placed in a group!");
end

if sum(groupSizes) ~= size(nonPaddedIdx,2)
   error("missing indices!") ;
end
end
