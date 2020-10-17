function [] = StabilizedExplicit()

addpath('../')
% 1. Define the problem

T = 1E-3;


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.015;
% eSize = 0.5;

model = createpde(1);

dx = 0.1; dy = 1;
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
for nSteps = NSteps
    dt = t/nSteps;
    [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps);
    
    
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeAnalyticalSolutionAndNorms( Nodes, Elements, t, CP, GPInfo, U);
    
    figure(99)
    subplot(2,2,1)
    loglog( NSteps(1:i), L2, '*-.')
    xlabel('nSteps')
    ylabel('L2');
    %         hold on
    
    subplot(2,2,2)
    loglog( NSteps(1:i), LInf, '*-.')
    xlabel('nSteps')
    ylabel('LInf');
    %         hold on
    
    subplot(2,2,3)
    loglog( NSteps(1:i), L2U, '*-.')
    xlabel('nSteps')
    ylabel('L2 u');
    %         hold on
    
    subplot(2,2,4)
    loglog( NSteps(1:i), LInfU, '*-.')
    xlabel('nSteps')
    ylabel('LInf u');
    %         hold on
    
    i = i+1;
    
    
 
end




function [L2, L2U, LInf, LInfU] = ComputeAnalyticalSolutionAndNorms(Nodes, Elements, t, CP, GPInfo, Xnum)

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, 3, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


Xa = real(Xa);

[L2, L2U, LInf, LInfU] = ComputeErrorNorms(Xnum, Xa, Nodes, Elements, GPInfo, CP);

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

