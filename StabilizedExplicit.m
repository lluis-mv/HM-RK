function [] = StabilizedExplicit()

dt = 1E-2;

eSize = 0.05;
model = createpde(1);

dx = 0.1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');


figure(1)
pdeplot(model)
drawnow


Nodes = mesh.Nodes';
Elements = mesh.Elements';

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);




[ElementMatrices] = ComputeElementalMatrices(Nodes, Elements);

[C, K ] = EnsambleMatrices(Nodes, Elements, ElementMatrices, dt);

[C, K, X0] = ApplyBC(Nodes, C, K);

A = C\(-K);
hola = 1;


ii = eye(3*nNodes, 3*nNodes);



B = ii+dt*A;

X = X0;


b = eig(B);
figure(900)
plot(sort(b));

for t = 1:50
    X = B*X;
end

hola = 1;


dofsWP = 3*([1:nNodes]-1)+3;
figure(1)
s = trimesh(Elements, Nodes(:,1), Nodes(:,2), X(dofsWP), 'FaceColor', 'interp', 'EdgeColor', 'none');
view(0, 90)
colorbar
axis equal
hola = 1;
function [C, K, X0] = ApplyBC(Nodes, C, K)

nNodes = size(Nodes, 1);

nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));

% Fix wp on top
dofs = 3*(nodesTop-1)+3;

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

% Fix uY bottom
dofs = 3*(nodesBottom-1)+2;

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

% Fix uX on left and Right
dofs = 3*([nodesLeft; nodesRight]-1)+1;

C(dofs,:) = 0;
K(dofs,:) = 0;
C(dofs,dofs) =eye(length(dofs));

X0 = zeros(3*nNodes, 1);
for i = 1:nNodes
    X0(3*(i-1)+3) = 1;
end
% Fix wp on top
dofs = 3*(nodesTop-1)+3;
X0(dofs) = 0;

function [C, K] = EnsambleMatrices(Nodes, Elements, ElementMatrices, dt)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = sparse(nSystem, nSystem);
K = C;

perme = 1E-3;
one = [1,1,0]';


De = zeros(6,6);

E = 1000;
nu = 0.3;

for i = 1:3
    for j = 1:3
        if ( i==j)
            De(i,j) = (1-nu);
        else
            De(i,j) = nu;
        end
    end
    De(i+3, i+3) = (1-2*nu)/2;
end

De = E/(1+nu)/(1-2*nu) * De;


De = De([1,2,4], [1,2,4]);
M=E*(1-nu)/(1+nu)/(1-2*nu);

for el = 1:nElements
    ind = Elements(el,:);
    index = [];
    for ii = 1:length(ind)
        index = [ index, (ind(ii)-1)*3 + [1,2,3] ];
    end
    
    kke = ElementMatrices(el).B'*De*ElementMatrices(el).B;
    Q = ElementMatrices(el).B'*one * ElementMatrices(el).N;
    H = -ElementMatrices(el).dN_dX'*perme*ElementMatrices(el).dN_dX;
    
    he = ElementMatrices(el).he;
    AlphaStab = 3/M+6*perme*dt/he^2;
    
    AlphaStab = 3/M + 6*dt*perme/he^2;
    AlphaStab = -AlphaStab*1.5;
    
    Ms = ElementMatrices(el).Ms * AlphaStab;
    Ce = [kke, Q; Q', Ms];
    Ke = [zeros(6,9); zeros(3,6), H];
    
    aux = [1,2,7,3,4,8,5,6,9];
    
    Ke = Ke(aux,aux);
    Ce = Ce(aux,aux);
    
    K(index,index) =  K(index,index) + Ke*ElementMatrices(el).Weight;
    C(index,index) =  C(index,index) + Ce*ElementMatrices(el).Weight;
    
end




function [ElementMatrices] = ComputeElementalMatrices(Nodes, Elements)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
ndim = 2;

for el = 1:nElements
    
    X = Nodes(Elements(el,:),:);
    % Linear triangles
    alfa = 1/3;
    beta = 1/3;
    
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
    Nsmall_chi = [-1 -1; 1 0; 0 1];
    Nu = (zeros(ndim, nNodes*ndim));
    for i = 1:3
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
    J = Nsmall_chi'*X;
    dN_dX = inv(J)*Nsmall_chi';
    
    B = [];
    for i = 1:3
        b = [-dN_dX(1,i), 0; 0, -dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
        B = [B, b];
    end
    
    Area = [1 1 1;
        X(1,1) X(2,1) X(3,1);
        X(1,2) X(2,2) X(3,2)];
    Area = det(Area)/2;
    
    he = 0;
    for i = 1:3
        aux = 0;
        for j = 1:2
            aux = aux + dN_dX(j,i);
        end
        he = he + abs(aux);
    end
    he = sqrt(2)*he;
    he = 4/he;
    Ms = 1/18*[2,-1,-1;-1,2,-1;-1,-1,2];
    ElementMatrices(el).Weight = Area;
    ElementMatrices(el).B =B;
    ElementMatrices(el).dN_dX = dN_dX;
    ElementMatrices(el).N = Nsmall';
    ElementMatrices(el).Nu = Nu;
    ElementMatrices(el).he = he;
    ElementMatrices(el).Ms = Ms;
end
