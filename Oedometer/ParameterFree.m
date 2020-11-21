function [] = OedometerElasticRK()

addpath('../')
% 1. Define the problem

T = 1E-5;


CP.E = 1000;
CP.nu = 0.3;
CP.k = 1E-2;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.02;

model = createpde(1);

dx = 0.05; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';



ElementType = 'T3T3';
Nodes = Nodes;
Elements = Elements;

NSteps = 10.^[0:5];
NSteps = floor(NSteps); NSteps = sort(NSteps);





figure(99)
for i = 1:4; subplot(2,2,i); hold off; end


RKMethod = 1;

kk = 10.^[-5:2:-1];
EE = 10.^[3:1:5];

for k = kk
    for E = EE
        
        i = 1;
        
        
        CP.E = E;
        nu = CP.nu;
        CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu)
        CP.k = k;
        
        t = T/CP.M/CP.k;
        
        for nSteps = NSteps
            
            dt = t/nSteps;
            
            [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod);
            
            
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
            
            
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP);
            
            
            i = i+1;
            
        end
        figure(99)
        subplot(2,2,1)
        loglog( NSteps, L2, '*-.')
        xlabel('nSteps')
        ylabel('L2');
        hold on
        
        
        subplot(2,2,2)
        loglog( NSteps, LInf, '*-.')
        xlabel('nSteps')
        ylabel('LInf');
        hold on
        
        
        subplot(2,2,3)
        loglog( NSteps, L2U, '*-.')
        xlabel('nSteps')
        ylabel('L2 u');
        hold on
        
        
        subplot(2,2,4)
        loglog( NSteps, LInfU, '*-.')
        xlabel('nSteps')
        ylabel('LInf u');
        hold on
    end
end




function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 3, 3, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


Xa = real(Xa);



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



% LInfU = LInfU*CP.M;
% L2U = L2U*CP.M;

