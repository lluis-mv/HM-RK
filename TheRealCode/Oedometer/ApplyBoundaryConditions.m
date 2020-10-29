

function [C, K, X0, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));

% Fix wp on top
dofs = 3*(nodesTop-1)+3;
nDirichlet = [nDirichlet; dofs];


% dofs = 3*( [1:nNodes ]-1)+3;
C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

% Fix uY bottom
dofs = 3*(nodesBottom-1)+2;
nDirichlet = [nDirichlet; dofs];

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

% Fix uX on left and Right
dofs = 3*([nodesLeft; nodesRight]-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

X0 = zeros(3*nNodes, 1);
for i = 1:nNodes
    X0(3*(i-1)+3) = 1;
end
% Fix wp on top
dofs = 3*(nodesTop-1)+3;
nDirichlet = [nDirichlet; dofs];
X0(dofs) = 0;

f = zeros(3*nNodes, 1);

for el = 1:nElements
    Cel = Elements(el,:);
    found = false;
    for i = 1:3
        for j = i+1:3
            if ( any(Cel(i) == nodesTop))
                if ( any(Cel(j) == nodesTop))
                    found = true;
                    ii = i;
                    jj = j;
                end
            end
        end
    end
    if (found)
        nodi = Cel(ii);
        nodj = Cel(jj);
        XX = Nodes(nodi,:)-Nodes(nodj,:);
        
        normal = [XX(2), -XX(1)];
        normal = normal/norm(normal);
        fe = 0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
        
        index = [ 3*(nodi-1)+[1,2], 3*(nodj-1)+[1,2]];
        f(index) = f(index) + fe;
        
    end
end

f = 0*f;
fini = f;
