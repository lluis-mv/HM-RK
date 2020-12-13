
function [C, K, X0, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K)



penalty = 1;

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == min(Nodes(:,2)));
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));



% Fix wp on top
index = find( (Nodes(:,1) > 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*(index-1)+3;
% dofs = 3*([1:nNodes]'-1)+3;
nDirichlet = [nDirichlet; dofs];

K(dofs,:) = 0;
C(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));



% Fix uY bottom
index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([nodesBottom; index]-1)+2;
nDirichlet = [nDirichlet; dofs];

K(dofs,:) = 0;
C(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));

% Fix uX on left and Right
dofs = 3*([index; nodesLeft; nodesRight]-1)+1;
nDirichlet = [nDirichlet; dofs];
K(dofs,:) = 0;
C(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));




X0 = zeros(3*nNodes, 1);


X0(dofs) = 0;







f = zeros(3*nNodes, 1);
indexL = find( (Nodes(:,1) >= 1) &  abs(Nodes(:,2)) < 1E-8);

if (size(Elements,2) == 3)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = false;
        for i = 1:3
            for j = i+1:3
                if ( any(Cel(i) == indexL))
                    if ( any(Cel(j) == indexL))
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
elseif (size(Elements,2) == 6)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = 0;
        indexs = [];
        for i = 1:6
            if ( any(Cel(i) == indexL))
                found = found+1;
                indexs = [indexs, Cel(i)];
            end
        end
        if (found == 3)
            nodi = (indexs(1));
            nodj = (indexs(2));
            XX = Nodes(nodi,:)-Nodes(nodj,:);
            
            
            normal = [XX(2), -XX(1)];
            normal = normal/norm(normal);
            ff = [ 1/6,   0, 1/6,   0, 2/3,   0;
                    0, 1/6,   0, 1/6,   0, 2/3]';
            fe = 1*ff*normal'*norm(XX);

            
            index = [];
            for i = 1:3
                index = [index, 3*(indexs(i)-1)+[1,2]];
            end
            f(index) = f(index) + fe;
        end
    end
    
end

fini = f;

fini = 0*f;

