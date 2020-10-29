function [] = HelpMe2()

% load('thisFile.mat')
% x0 = x;

x0 = [9];


x0 = [9];
ComputeThisErrors(x0)
hold on

x0 = [7];
ComputeThisErrors(x0)
hold on

x0 = [5];
ComputeThisErrors(x0)
hold on


x0 = [4.5];
ComputeThisErrors(x0)
hold off

return

fun = @(x) ComputeThisErrors(x);
options = optimset('PlotFcns',@optimplotfval);
x = fminsearch(fun,x0,options)
save('thisFile', 'x')


function [L2] = ComputeThisErrors(x)
x

a = 0;
b = x(1);


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-1;
CP.a = a;
CP.b = b;

i = 1;
Ni = [1,10,100,1000];
for n = Ni
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeThisProblem(n, CP,1E-1,1);
    i = i+1;
end
for n = Ni
    CP.E = 1000;
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeThisProblem(n, CP,1E-3,10);
    i = i+1;
end
N = [Ni, 10*Ni(end)*Ni];

for n = Ni
    CP.k = 1E-4;
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeThisProblem(n, CP,1E-5,5);
    i = i+1;
end
for n = Ni
    CP.k = 1E-6;
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeThisProblem(n, CP,1E-7,0.5);
    i = i+1;
end
N = [N, 10*N(end)*N];







figure(900)
loglog(N, L2, 'b*-.',N, LInf, 'c*-.',N, L2U, 'r*-.',N, LInfU, 'm*-.')
xlabel('nSteps')
ylabel('error norms')




index = find(isnan(L2));

ind2 = find(L2 < 10E500);
L2 = sum(L2(ind2))*10^length(index);







function [L2, L2U, LInf, LInfU] = ComputeThisProblem(nSteps, CP, T, H)

if (nargin == 2)
    T = 1E-1;
end
% nSteps = 10;

if (nargin == 1)
    CP.E = 100;
    CP.nu = 0.3;
    CP.k = 1E-3;
end

nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;
dt = T/nSteps/CP.M/CP.k;

eSize = H/30;
model = createpde(1);

dx = H/10; 
dy = H;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');


% figure(1)
% pdeplot(model)
% drawnow


Nodes = mesh.Nodes';
Elements = mesh.Elements';

nNodes = size(Nodes, 1)
nElements = size(Elements, 1);




[ElementMatrices] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, ElementMatrices, CP, dt);

[C, K, X0, f] = ApplyBC(Nodes, Elements, C, K);

A = C\(K);
hola = 1;

ii = eye(3*nNodes, 3*nNodes);
B = ii+dt*A;
% Second order RK
B = ii+0.5*dt*A+0.5*dt*A*(ii+dt*A);
b = eig(B);

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, ElementMatrices, CP, dt, 0);

[Ca, Ka, ~, ~] = ApplyBC(Nodes, Elements, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

XA = 0*X0;

c = (vectors)\X0;

for i = 1:size(values, 1)
    XA = XA + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


XA = real(XA);


ii = eye(3*nNodes, 3*nNodes);
B2 = ii+dt*Aa;
b2 = eig(B2);

% figure(900)
% plot(sort(b2), 'r')
% hold on
% plot( sort(b), 'k-.');
% hold off;
% ylim([-1.5,1.5])




f = 0*f; X= X0;





for i = 1:nSteps
    X = B*X;
end


[L2, L2U, LInf, LInfU] = ComputeErrorNorms(X, XA, Nodes, Elements, ElementMatrices, CP);


% %  implicit integration
% B = ii-dt*A;
% X = zeros(size(X));
% X = B\(0*X+ C\fa);




dofsWP = 3*([1:nNodes]-1)+3;
% figure(1)
% s = trimesh(Elements, Nodes(:,1), Nodes(:,2), X(dofsWP), 'FaceColor', 'interp', 'EdgeColor', 'none');
% view(0, 90)
% colorbar
% axis equal
% hola = 1;
% 
% 
% figure(2)
% s = trimesh(Elements, Nodes(:,1), Nodes(:,2), X(dofsWP)-XA(dofsWP), 'FaceColor', 'interp', 'EdgeColor', 'none');
% view(0, 90)
% colorbar
% axis equal
% hola = 1;

function [L2, L2U, LInf, LInfU] = ComputeErrorNorms(X, Xa, Nodes, Elements, ElementMatrices, CP)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


indexWP = 3*[1:nNodes];
LInf = max( abs( X(indexWP)-Xa(indexWP)));
LInfU = 0;
for i = 1:nNodes
    ind = 3*(i-1)+[1,2];
    thisNorm = norm(Xa(ind)-X(ind));
    LInfU = max(LInfU, thisNorm);
end

L2 = 0;
L2U = 0;


alfa = 2/3; beta = 1/6;
N1 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 1/6;
N2 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 2/3;
N3 = [ 1 - alfa - beta, alfa,  beta];

for el = 1:nElements
    Cel = Elements(el,:);
    indWP = 3*(Cel-1)+3;
    err = ( Xa(indWP)-X(indWP));
    L2 = L2 + ElementMatrices(el).Weight/3* ( abs(N1*err) + abs(N2*err)+abs(N3*err));
    
    indx = 3*(Cel-1)+1;
    indy = 3*(Cel-1)+2;
    ux = Xa(indx)-X(indx);
    uy = Xa(indy)-X(indy);
    L2U = L2U + ElementMatrices(el).Weight/3* ( norm(N1*[ux,uy]) + norm(N2*[ux,uy])+norm(N3*[ux,uy]));
    
end



LInfU = LInfU*CP.M;
L2U = L2U*CP.M;




function [C, K, X0, f] = ApplyBC(Nodes, Elements, C, K)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));

% Fix wp on top
dofs = 3*(nodesTop-1)+3;

% dofs = 3*( [1:nNodes ]-1)+3;

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
dofs = 3*([1:nNodes]-1)+1;
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


function [C, K] = EnsambleMatrices(Nodes, Elements, ElementMatrices, CP, dt, AlphaStabM)

if (nargin == 5)
    AlphaStabM = 1;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = zeros(nSystem, nSystem);
K = C;


one = [1,1,0]';

perme = CP.k;
ConstModulus=CP.M;

for el = 1:nElements
    ind = Elements(el,:);
    index = [];
    for ii = 1:length(ind)
        index = [ index, (ind(ii)-1)*3 + [1,2,3] ];
    end
    
    kke = ElementMatrices(el).B'*ElementMatrices(el).D*ElementMatrices(el).B;
    Q = ElementMatrices(el).B'*one * ElementMatrices(el).N;
    H = -ElementMatrices(el).dN_dX'*perme*ElementMatrices(el).dN_dX;
    
    he = ElementMatrices(el).he;
    he = sqrt(ElementMatrices(el).Weight);
    AlphaStab = CP.a/ConstModulus+CP.b*perme*dt/he^2;
    
    if ( AlphaStab < 0)
        AlphaStab = 0;
    end
    
    AlphaStab = AlphaStab*AlphaStabM;
    
    
    Ms = ElementMatrices(el).Ms * AlphaStab;
    Ce = [kke, Q; -Q', Ms];
    Ke = [zeros(6,9); zeros(3,6), H];
    
    aux = [1,2,7,3,4,8,5,6,9];
    
    Ke = Ke(aux,aux);
    Ce = Ce(aux,aux);
    
    K(index,index) =  K(index,index) + Ke*ElementMatrices(el).Weight;
    C(index,index) =  C(index,index) + Ce*ElementMatrices(el).Weight;
    
end




function [ElementMatrices] = ComputeElementalMatrices(Nodes, Elements, CP)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
ndim = 2;

De = zeros(6,6);

E = CP.E;
nu = CP.nu;

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
    ElementMatrices(el).D = De;
end
