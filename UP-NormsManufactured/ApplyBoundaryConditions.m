
function [C, K, X0, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K)



penalty = 1;

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));

% Fix wp on top
% dofs = 3*([1:nNodes]'-1)+3;
dofs = [];



nDirichlet = [nDirichlet; dofs];


% dofs = 3*( [1:nNodes ]-1)+3;
C(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

% Fix uY bottom
dofs = 3*([nodesTop; nodesBottom]-1)+2;
% dofs = 3*([1:nNodes]'-1)+2;
nDirichlet = [nDirichlet; dofs];

C(dofs,:) = 0;

C(dofs,dofs) = penalty*eye(length(dofs));

% Fix uX on left and Right

dofs = 3*( [nodesLeft; nodesRight]-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;

C(dofs,dofs) =penalty*eye(length(dofs));

X0 = zeros(3*nNodes, 1);


% Fix wp on top
dofs = 3*(nodesTop-1)+3;
X0(dofs) = 0;



fini = 0*X0;








