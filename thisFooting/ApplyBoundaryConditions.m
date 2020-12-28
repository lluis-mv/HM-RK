
function [C, K, X0, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K)



penalty = 1;

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == min(Nodes(:,2)));
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));



dofs = 3*([1:nNodes]'-1)+3;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));


% Fix uY bottom
index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([nodesBottom; index]-1)+2;
nDirichlet = [nDirichlet; dofs];

C(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));

% Fix uX on left and Right

dofs = 3*([nodesBottom; nodesRight; nodesLeft]-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

X0 = zeros(3*nNodes, 1);


fini = 0*X0;




