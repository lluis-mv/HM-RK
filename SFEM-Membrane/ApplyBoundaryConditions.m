
function [C, K, X0, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K)

if (nargin ~= 5)
    error('it should be five!!!!')
end

penalty = 1;

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nDirichlet = [];

nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));

% Fix wp on top
dofs = 3*(nodesLeft-1)+3;
nDirichlet = [nDirichlet; dofs];
% dofs = 3*([1:nNodes]'-1)+3;
% nDirichlet = [nDirichlet; dofs];



C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

% Fix uY bottom
dofs = 3*(nodesLeft-1)+2;
nDirichlet = [nDirichlet; dofs];

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));

% Fix uX on left and Right
dofs = 3*([nodesLeft; ]-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

X0 = zeros(3*nNodes, 1);




% Fix wp on top
dofs = 3*(nodesLeft-1)+3;
X0(dofs) = 0;


% now try that
if ( length([GPInfo(1,1).dofsWP]) ~= length([GPInfo(1,1).dofsWPreal]) )
    for el = 1:nElements
        dofsWP = GPInfo(el,1).dofsWP;
        dofsReal = GPInfo(el,1).dofsWPreal;
        
        KK = 1/2*[1,1,0;
            0, 1, 1;
            1, 0,1];
        X0( dofsReal(4:6)) = KK*X0(dofsWP) ;
          
    end
end


f = zeros(3*nNodes, 1);


if (size(Elements,2) == 3)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = false;
        for i = 1:3
            for j = i+1:3
                if ( any(Cel(i) == nodesRight))
                    if ( any(Cel(j) == nodesRight))
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
            
            normal = [XX(1), XX(2)];
            normal = normal/norm(normal);
            fe = 0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
            
            index = [ 3*(nodi-1)+[1,2], 3*(nodj-1)+[1,2]];
            f(index) = f(index) + fe;
            fe
        end
    end
elseif (size(Elements,2) == 6)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = 0;
        indexs = [];
        for i = 1:6
            if ( any(Cel(i) == nodesRight))
                found = found+1;
                indexs = [indexs, Cel(i)];
            end
        end
        if (found == 3)
            nodi = (indexs(1));
            nodj = (indexs(2));
            XX = Nodes(nodi,:)-Nodes(nodj,:);
            
            
            normal = [XX(1), XX(2)];
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
f = 0*f;
