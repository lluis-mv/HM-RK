
function [C, K, X0, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K)

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
dofs = 3*([nodesTop; nodesBottom]-1)+3;



nDirichlet = [nDirichlet; dofs];


% dofs = 3*( [1:nNodes ]-1)+3;
C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

% Fix uY bottom
dofs = 3*([nodesTop; nodesBottom]-1)+2;
% dofs = 3*([1:nNodes]'-1)+2;
nDirichlet = [nDirichlet; dofs];

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) = penalty*eye(length(dofs));

% Fix uX on left and Right
dofs = 3*([1:nNodes]'-1)+1;
nDirichlet = [nDirichlet; dofs];
C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =penalty*eye(length(dofs));

X0 = zeros(3*nNodes, 1);
for i = 1:nNodes
    X0(3*(i-1)+3) = 1;
end

K(nDirichlet, nDirichlet) = eye(length(nDirichlet));


% Fix wp on top
dofs = 3*(nodesTop-1)+3;
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





for el = 1:nElements
    Cel = Elements(el,:);
    dofsU = GPInfo(el,1).dofsU;
    dofsWP = GPInfo(el,1).dofsWP;
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:size(GPInfo, 2)
        
        
        
        Xpg = GPInfo(el, gp).Nu(1,1:2:end)*Nodes(Cel,:);
        y = Xpg(2);
        ff =[0 ; (y*(2*y - 2))/25 + (y - 1)^2/50 + y^2/50];
        ff2 = (y*(y - 1)^2)/50 + (y^2*(2*y - 2))/100 - 2;
        
        f( GPInfo(el, 1).dofsU) = f( GPInfo(el, 1).dofsU) - GPInfo(el, gp).Weight*GPInfo(el, gp).Nu'*( ff);
        
        f( GPInfo(el, 1).dofsWP) = f( GPInfo(el, 1).dofsWP) - GPInfo(el, gp).Weight*GPInfo(el, gp).N'*ff2;
    end
    
end


fini = 0*f;

X0 = 0*X0;
hola = 1;





