
function [] = OedometerElastic()

addpath('../')
% 1. Define the problem

T = 4E-3;


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-3;

nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.015;
eSize = 0.015;
eSize = 0.03;
eSize = 0.04;
% eSize = 0.0075;
model = createpde(1);

dx = 0.1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';




NSteps = 10.^linspace(0, 5, 10);
NSteps = 10.^[0:4];
NSteps = floor(NSteps); NSteps = sort(NSteps);

i = 1;

firstTime = true;
% NSteps = 1000;

L2 = NSteps*nan; L2U = NSteps*nan; LInf = NSteps*nan; LInfU = NSteps*nan;
L2i = NSteps*nan; L2Ui = NSteps*nan; LInfi = NSteps*nan; LInfUi = NSteps*nan;
L2g = NSteps*nan; L2Ug = NSteps*nan; LInfg = NSteps*nan; LInfUg = NSteps*nan;

for nSteps = NSteps
    dt = t/nSteps;
    
    [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', 3);
    figure(904)
    plot(U(3:3:end), Nodes(:,2), 'b*')
    hold on
    
    figure(905)
    plot(U(2:3:end), Nodes(:,2), 'b*')
    hold on
    
    
    if ( firstTime)
        [Xa] = ComputeAnalyticalSolution(Nodes, Elements, 'T3T3', t, CP, GPInfo,U);
    end
    figure(905)
    hold on
    plot( Xa(2:3:end), Nodes(:,2), 'bs')
    figure(904)
    hold on
    plot( Xa(3:3:end), Nodes(:,2), 'bs')
    drawnow
    [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP);
    
    
    
    [U,GPInfo] = ComputeThisLinearProblem(Nodes2, Elements2, CP, dt, nSteps, 'T6T6', 1);
    if ( firstTime)
        [Xa2] = ComputeAnalyticalSolution(Nodes2, Elements2,'T6T6',  t, CP, GPInfo,U);
    end
    
    figure(905)
    hold on
    plot( Xa2(2:3:end), Nodes2(:,2), 'rs')
    figure(904)
    hold on
    plot( Xa2(3:3:end), Nodes2(:,2), 'rs')
    [L2i(i), L2Ui(i), LInfi(i), LInfUi(i)] = ComputeErrorNorms(U, Xa2, Nodes2, Elements2, GPInfo, CP);
    
    figure(904)
    hold on
    plot(U(3:3:end), Nodes2(:,2), 'r*')
    
    
    figure(905)
    hold on
    plot(U(2:3:end), Nodes2(:,2), 'r*')
    
    
    
    [U,GPInfo] = ComputeThisLinearProblem(Nodes2, Elements2, CP, dt, nSteps, 'T6T3', 3);
    if ( firstTime)
        [Xa3] = ComputeAnalyticalSolution(Nodes2, Elements2,'T6T3',  t, CP, GPInfo,U);
        firstTime = false;
    end
    [L2g(i), L2Ug(i), LInfg(i), LInfUg(i)] = ComputeErrorNorms(U, Xa3, Nodes2, Elements2, GPInfo, CP);
    
    figure(905)
    hold on
    plot( Xa3(2:3:end), Nodes2(:,2), 'gs')
    figure(904)
    hold on
    plot( Xa3(3:3:end), Nodes2(:,2), 'gs')
    
    figure(904)
%     hold on
%     plot(U(3:3:end), Nodes2(:,2), 'g*')
    hold off
%     
    figure(905)
%     hold on
%     plot(U(2:3:end), Nodes2(:,2), 'g*')
    hold off
%     
    
    
    
    figure(99)
    subplot(2,2,1)
    loglog( NSteps, L2, 'b*-.', NSteps, L2i, 'r*-.',  NSteps, L2g, 'g*-.')
    xlabel('nSteps')
    ylabel('L2');
    
    
    subplot(2,2,2)
    loglog( NSteps, LInf, 'b*-.', NSteps, LInfi, 'r*-.', NSteps, LInfg, 'g*-.')
    xlabel('nSteps')
    ylabel('LInf');
    
    
    subplot(2,2,3)
    loglog( NSteps, L2U, 'b*-.', NSteps, L2Ui, 'r*-.', NSteps, L2Ug, 'g*-.')
    xlabel('nSteps')
    ylabel('L2 u');
    
    
    subplot(2,2,4)
    loglog( NSteps, LInfU, 'b*-.',  NSteps, LInfUi, 'r*-.',  NSteps, LInfUg, 'g*-.')
    xlabel('nSteps')
    ylabel('LInf u');
    
    
    i = i+1;
    
end




function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
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




for nod = 1:nNodes
    z= 1-Nodes(nod,2);
    
    uu = z-1;
    
    for m = 0:100
        
        term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
        uu = uu+term;
    end
    Xa(3*(nod-1)+2) = uu/M;
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

return;

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

