
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

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;

wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;

auxK = [Na2, Na2, wa;
    Na1, Na2, wa;
    Na2, Na1, wa;
    Nb2, Nb2, wb ;
    Nb1, Nb2, wb ;
    Nb2, Nb1, wb ];

al = auxK(:,1)';
be = auxK(:,2)';
w = auxK(:,3)'/sum(auxK(:,3)');

f = zeros(3*nNodes,1);
for el = 1:nElements
    Cel = Elements(el,:);
    dofsU = GPInfo(el,1).dofsU;
    dofsWP = GPInfo(el,1).dofsWPreal;
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:length(w)
        [Nu, Np] = GetShapeFunctions( al(gp), be(gp), length(dofsU), length(dofsWP));
        
        
        Xpg = Nu(1,1:2:end)*Nodes(Cel,:);
         y = Xpg(2);
         ff2 =1/50;
         ff = y/50-1/100;
        f( GPInfo(el, 1).dofsWP) = f( GPInfo(el, 1).dofsWP) - wA*w(gp)*Np'*( ff);
        
        f( GPInfo(el, 1).dofsWP-1) = f( GPInfo(el, 1).dofsWP-1) - wA*w(gp)*Np'*( ff2);
    end
    
end


fini = 0*f;

X0 = 0*X0;
hola = 1;









function [Nu, Np] = GetShapeFunctions( alfa, beta, nU, nP)

ndim = 2;
if (nU == 6)
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
    
    Nu = (zeros(ndim, 3*ndim));
    for i = 1:3
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
elseif (nU == 12)
    Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];
    
    
    Nu = (zeros(ndim, 6*ndim));
    for i = 1:6
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
end


if ( nP == 3)
    Np =  [ 1 - alfa - beta; alfa;  beta]';
elseif ( nP == 6)
    Np =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)]';
end




