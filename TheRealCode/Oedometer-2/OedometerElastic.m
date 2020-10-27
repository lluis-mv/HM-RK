function [] = OedometerElastic()

addpath('../')
% 1. Define the problem

T = 1E-8;


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.015;
eSize = 0.02;
% eSize = 0.0075;
model = createpde(1);

dx = 0.01; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

% figure(1)
% pdeplot(model)
% drawnow


Nodes = mesh.Nodes';
Elements = mesh.Elements';





NSteps = 10.^linspace(0, 5, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);

i = 1;

firstTime = true;

for nSteps = NSteps
    dt = t/nSteps;
    
    [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps);   
    
    if ( firstTime)
        [Xa] = ComputeAnalyticalSolution(Nodes, Elements, t, CP, GPInfo,U);
        firstTime = false;
    end
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP);
    
    [U,GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps);   
    [L2i(i), L2Ui(i), LInfi(i), LInfUi(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP);
    
    
    figure(99)
    subplot(2,2,1)
    loglog( NSteps(1:i), L2, 'b*-.', NSteps(1:i), L2i, 'r*-.')
    xlabel('nSteps')
    ylabel('L2');
    
    
    subplot(2,2,2)
    loglog( NSteps(1:i), LInf, 'b*-.', NSteps(1:i), LInfi, 'r*-.')
    xlabel('nSteps')
    ylabel('LInf');
    
    
    subplot(2,2,3)
    loglog( NSteps(1:i), L2U, 'b*-.', NSteps(1:i), L2Ui, 'r*-.')
    xlabel('nSteps')
    ylabel('L2 u');
    
    
    subplot(2,2,4)
    loglog( NSteps(1:i), LInfU, 'b*-.',  NSteps(1:i), LInfUi, 'r*-.')
    xlabel('nSteps')
    ylabel('LInf u');
    
    
    i = i+1;
    
end




function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;
% Other analytical solution...
nNodes = size(Nodes,1);
M = CP.M;
k = CP.k;
for nod = 1:nNodes
    xx = 1-Nodes(nod,2);
    TT = M * t*k;
    pw = 0;
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw = pw + 2/aux * sin( aux * xx) * exp( - aux^2 * TT);
    end
    Xa(3*(nod-1)+3) = pw;
end




function [L2, L2U, LInf, LInfU] = ComputeErrorNorms(X, Xa, Nodes, Elements, GPInfo, CP)



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
    L2 = L2 + GPInfo(el).Weight/3* ( abs(N1*err) + abs(N2*err)+abs(N3*err));
    
    indx = 3*(Cel-1)+1;
    indy = 3*(Cel-1)+2;
    ux = Xa(indx)-X(indx);
    uy = Xa(indy)-X(indy);
    L2U = L2U + GPInfo(el).Weight/3* ( norm(N1*[ux,uy]) + norm(N2*[ux,uy])+norm(N3*[ux,uy]));
    
end



LInfU = LInfU*CP.M;
L2U = L2U*CP.M;

