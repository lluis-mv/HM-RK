
function [C, K, X0, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K, dofsPerNode)

if (nargin == 5)
    dofsPerNode = 3;
end
if (dofsPerNode~=3)
    error('You have to code it')
end


penalty = 1;

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == min(Nodes(:,2)));
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));


index = find( (Nodes(:,1) > 1) &  abs(Nodes(:,2)) < 1E-8);
index = [nodesBottom];

% ho he de fer millor per separar allò....
% i.e. the element is T6T3
if ( length([GPInfo(1,1).dofsWP]) ~= length([GPInfo(1,1).dofsWPreal]) )
    dofsR = [];
    for el = 1:nElements
        nF = 0;
        ind = [];
        for ij = 1:6
            if ( any(Elements(el,ij) == index) )
                nF = nF+1;
                ind = [ind, Elements(el,ij)];
            end
        end
        if ( nF == 3)
            dofsR = [dofsR, ind];
        end
    end
    dofsR = unique(dofsR)';
    dofs = 3*(dofsR-1)+3;
else
    dofs = 3*([index]-1)+3;
end



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

dofs = 3*([nodesBottom; index]-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

X0 = zeros(3*nNodes, 1);
% X0(3:3:end) = -10;



fini = 0*X0;




