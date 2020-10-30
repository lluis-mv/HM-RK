function []=OneDimensionalEqual()
addpath('../ModifiedCamClay/')

nodes = linspace(0, 1, 25);


% dt = 0.01;
% nSteps = 200;


for i = 1:4
    figure(i);
    hold off;
end

M = 1;
k = 1;
for RKMethod = [ 1,2,3,4,5,6,7,8]
    
    
    TT = 1E-4;
    NSteps = 10.^linspace(0, 5, 10);
    NSteps = floor(NSteps);
    

    i = 1;
    for nSteps = NSteps
        dt = TT/nSteps;
        [L2(i), LInf(i), L2DISPL(i), LInfDISPL(i)] = CalculateProblemAndNorms(nodes, nSteps, dt, RKMethod, M, k);
        
        
%         figure(1)
%         loglog( NSteps(1:i), L2, '*-.')
%         xlabel('nSteps')
%         ylabel('L2');
%         hold on
        
        figure(2)
        loglog( NSteps(1:i), LInf, '*-.')
        xlabel('nSteps')
        ylabel('LInf');
        hold on
        
%         figure(3)
%         loglog( NSteps(1:i), L2DISPL, '*-.')
%         xlabel('nSteps')
%         ylabel('L2 u');
%         hold on
        
        figure(4)
        loglog( NSteps(1:i), LInfDISPL, '*-.')
        xlabel('nSteps')
        ylabel('LInf u');
        hold on
        i = i+1;
    end
    
    clear L2
    clear LInf
    clear L2DISPL
    clear LInfDISPL
end




function [L2, LInf, L2DISPL, LInfDISPL ] = CalculateProblemAndNormsImplicit(nodes, nSteps, dt, RKMethod, M, k)

% M = 1;


% M = 1e-3;
% k = 1;
dt = dt/(M*k);


substepping = false;
if ( RKMethod > 10)
    alfa = floor(RKMethod/10);
    substepping = true;
end


h = nodes(2)-nodes(1);



AlphaStab = 0;
if ( dt < h^2/M/k/6)
    AlphaStab = 2/M-12*k*dt/h^2;
end


[AMatrix, FFExt, X0] = ConstructMatrices(nodes, M, k, AlphaStab, dt);
B = inv(eye(size(AMatrix))-dt*AMatrix);
b = eig(B);
min(b)
max(b)
figure(900)
plot(sort((b)), 'r')
[AMatrixA] = ConstructMatrices(nodes, M, k, 0, dt);
figure(900)
B = inv(eye(size(AMatrix))-dt*AMatrixA);
hold on
b = eig(B);
plot(sort((b)), 'k')
ylim([-1.2, 1.2])
drawnow
hold off
drawnow
[X, NodalWaterPressure, NodalDisplacement] = IntegrateProblemImplicit(AMatrix, FFExt, X0, RKMethod, dt, nSteps);

[NWA, NDA] = IntegrateAnalytical(X0, AMatrixA, dt*nSteps);

[L2, LInf, L2DISPL, LInfDISPL ] = ComputeNorms( nodes, M*k*dt*nSteps, NodalWaterPressure, NodalDisplacement, NWA, NDA, M);

L2DISPL = L2DISPL*M;
LInfDISPL = LInfDISPL*M;

figure(101)
hold on
plot(nodes, NWA, 'k:')
hold off


figure(102)
hold on
plot(nodes, NDA, 'k:')
hold off





function [L2, LInf, L2DISPL, LInfDISPL ] = CalculateProblemAndNorms(nodes, nSteps, dt, RKMethod, M, k)

% M = 1;


% M = 1e-3;
% k = 1;
dt = dt/(M*k);


substepping = false;
if ( RKMethod > 10)
    alfa = floor(RKMethod/10);
    substepping = true;
end


h = nodes(2)-nodes(1);


AlphaStab = 12*dt*k/h^2*(1-exp(- (1400*dt*k/h^2)^(2*RKMethod) ));




[AMatrix, FFExt, X0] = ConstructMatrices(nodes, M, k, AlphaStab, dt);
b = eig(eye(size(AMatrix))+dt*AMatrix);
min(b)
max(b)
figure(900)
plot(sort((b)), 'r')
[AMatrixA] = ConstructMatrices(nodes, M, k, 0, dt);
figure(900)
hold on
b = eig(eye(size(AMatrixA))+dt*AMatrixA);
plot(sort((b)), 'k')
ylim([-1.2, 1.2])
drawnow
hold off
drawnow
[X, NodalWaterPressure, NodalDisplacement] = IntegrateProblem(AMatrix, FFExt, X0, RKMethod, dt, nSteps);

[NWA, NDA] = IntegrateAnalytical(X0, AMatrixA, dt*nSteps);

[L2, LInf, L2DISPL, LInfDISPL ] = ComputeNorms( nodes, M*k*dt*nSteps, NodalWaterPressure, NodalDisplacement, NWA, NDA, M);

L2DISPL = L2DISPL*M;
LInfDISPL = LInfDISPL*M;

figure(101)
hold on
plot(nodes, NWA, 'k:')
hold off
ylabel('water pressure')


figure(102)
hold on
plot(nodes, NDA, 'k:')
hold off
ylabel('displacement')

function [NodalWaterPressure, NodalDisplacement] = IntegrateAnalytical(X0, AMatrix, t)

[vectors, values] = eig(full(AMatrix), 'nobalance');
%[vectors, values] = eig(full(AMatrix));
x = 0*X0;

c = (vectors)\X0;

for i = 1:size(values, 1)
    x = x + c(i)*exp(values(i,i)*t)*vectors(:,i);
end
hola = 1;

max(abs(imag(x)));
x = real(x);

nNodes = length(X0)/2;
NodalWaterPressure = (zeros(nNodes, 1));
NodalDisplacement = (zeros(nNodes, 1));
for i = 1:nNodes
    NodalWaterPressure(i) = x(2*(i-1)+2);
    NodalDisplacement(i) = x(2*(i-1)+1);
    
end

function [X, NodalWaterPressure, NodalDisplacement] = IntegrateProblem( AMatrix, FFExt, X, RKMethod, dt, nSteps)

[aRK, bRK] = GetRungeKutta(RKMethod);

k = zeros( length(X), length(bRK));

for t = 1:nSteps
    for i = 1:length(bRK)
        Xstep = X;
        for j = 1:i-1
            Xstep = Xstep + dt * aRK(i,j) * k(:,j);
        end
        
        [k(:,i)] = (AMatrix*Xstep)+FFExt;
        
    end
    
    Xnew = X;
    for i = 1:length(bRK)
        Xnew = Xnew + dt*bRK(i)*k(:,i);
    end
    X = Xnew;
    
    FFExt = 0*FFExt;
end

nNodes = length(X)/2;
NodalWaterPressure = (zeros(nNodes, 1));
NodalDisplacement = (zeros(nNodes, 1));
for i = 1:nNodes
    NodalWaterPressure(i) = X(2*(i-1)+2);
    NodalDisplacement(i) = X(2*(i-1)+1);
    
end


function [X, NodalWaterPressure, NodalDisplacement] = IntegrateProblemImplicit( AMatrix, FFExt, X, RKMethod, dt, nSteps)

one = eye(size(AMatrix));
A = one - dt*AMatrix;

invA = inv(A);
for t = 1:nSteps
   X = invA*X;
end

nNodes = length(X)/2;
NodalWaterPressure = (zeros(nNodes, 1));
NodalDisplacement = (zeros(nNodes, 1));
for i = 1:nNodes
    NodalWaterPressure(i) = X(2*(i-1)+2);
    NodalDisplacement(i) = X(2*(i-1)+1);
    
end



function [L2, LInf, L2DISPL, LInfDISPL ] = ComputeNorms( XX, TT, NodalWaterPressure, NodalDisplacement, NWA, NDA, M)

nGP = 3;
[xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
pwGP = interp1(XX, NWA, xGP);
L2 = EvaluateL2( XX, pwGP, NodalWaterPressure, nGP);

pw = EvalConsolidation(XX, TT);
LInf = max( abs(NodalWaterPressure - pw'));
LInf = max( abs(NodalWaterPressure - NWA));



[xGP, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP);
uuGP = interp1(XX, NDA, xGP);

L2DISPL = EvaluateL2( XX, uuGP, NodalDisplacement, nGP);

uu = EvalConsolidationDISPL(XX, TT);
uu = uu/M;
LInfDISPL = max( abs(NodalDisplacement - uu'));
LInfDISPL = max( abs(NodalDisplacement - NDA));

hola = 1;


figure(101)
plot(XX, pw, 'g', XX, NodalWaterPressure, 'r*')

figure(102)
plot(XX, uu, 'g', XX, NodalDisplacement, 'r*')

function [AMatrix, FFExt, X0] = ConstructMatrices(nodes, M, k, AlphaStab, dt)

nNodes = length(nodes);
CC = sparse(nNodes*2, nNodes*2);
KK = CC;
DeltaPW = -1;

for i = 1:nNodes-1
    h = nodes(i+1)-nodes(i);
    % ComputeElementalMatrix
    [Ke, Ce] = ComputeElemental(h, AlphaStab, M, k);
    
    % Ensamble Elemental Matrix
    CC( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        CC( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + Ce;
    KK( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        KK( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + Ke;
end

X0 = zeros(2*nNodes, 1);
for i = 1:nNodes
    X0(2*i) = 1;
end


% Apply dirichlet conditions
% ZeroWaterPressure
CC(2, :) = 0;
CC(2,2) = 1;
KK(2,:) = 0;



% Displacement
nn = 2*(nNodes -1)+1;
CC(nn, :) = 0;
CC(nn,nn) = 1;


%Right hand side (water pressure at dirichlet)
FFExt = (zeros(2*nNodes, 1));
FFExt(2) = DeltaPW/dt;

AMatrix = CC\KK;


FFExt = 0*FFExt;
X0(2) = 0;


function [xx, pw] = EvaluateSolution(dtEval, MEval)

% plot analytical solution
xx = linspace(0,1,100);
pw = 0*xx;
TT = MEval * dtEval;
for i = 1:length(xx)
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw(i) = pw(i) + 2/aux * sin( aux * xx(i)) * exp( - aux^2 * TT);
    end
end

function [Ke, Ce] = ComputeElemental(h, AlphaStab, M, k)

Ce =[  M/h, -M/h,               1/2,               1/2;
    -M/h,  M/h,              -1/2,              -1/2;
    -1/2,  1/2,  (AlphaStab*h)/12, -(AlphaStab*h)/12;
    -1/2,  1/2, -(AlphaStab*h)/12,  (AlphaStab*h)/12];

Ke = [ 0, 0,    0,    0;
    0, 0,    0,    0;
    0, 0,  k/h, -k/h;
    0, 0, -k/h,  k/h];

ind = [1,3,2,4];
Ce = Ce(ind,ind);
Ke = -Ke(ind,ind);




%%%% functions to evalute errors and so on...


function pw = EvalConsolidation(XX, TT)
pw = 0*XX;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:100
            aux = pi/2*(2*m+1);
            term = 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT);
            pw(i,j) = pw(i,j)+term ;
        end
    end
end

function uu = EvalConsolidationDISPL(XX, TT)
uu = -(1-XX);
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:100
            z=XX(i,j);
            term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
            uu(i,j) = uu(i,j)+term;
        end
    end
end
uu = -uu;
function [xPG, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP)

if ( nGP == 1)
    xPG = [];
    for i = 1:length(XX)-1
        xPG = [xPG; mean(XX(i:i+1))];
    end
elseif (nGP == 2)
    xPG=[];
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
elseif (nGP == 3)
    xPG=[];
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
end

uuGP = EvalConsolidationDISPL(xPG, TT);

function [xPG, pwPG] = EvalConsolidationGP(XX, TT, nGP)

if ( nGP == 1)
    xPG = [];
    for i = 1:length(XX)-1
        xPG = [xPG; mean(XX(i:i+1))];
    end
elseif (nGP == 2)
    xPG=[];
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
elseif (nGP == 3)
    xPG=[];
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
end


pwPG = EvalConsolidation(xPG, TT);

function subError = EvaluateL2( XX, pwGP, WPNodes, nGP)

subError = 0;

if ( nGP == 1)
    L = XX(2)-XX(1);
    w = L/2*2;
    for i = 1:length(XX)-1
        N = [0.5, 0.5];
        nodal = N*(WPNodes(i:i+1));
        subError = subError + (nodal-pwGP(i))^2*w;
    end
elseif (nGP == 2)
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    L = XX(2)-XX(1);
    w = L/2;
    for i = 1:length(XX)-1
        nodal = N * [WPNodes(i:i+1)];
        pwLocal = pwGP(i,:);
        for a = 1:2
            subError = subError + (nodal(a) -pwLocal(a))^2 * w;
        end
    end
elseif (nGP == 3)
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    L = XX(2)-XX(1);
    w = L/2*[5/9, 8/9,5/9];
    for i = 1:length(XX)-1
        nodal = N * [WPNodes(i:i+1)];
        pwLocal = pwGP(i,:);
        for a = 1:3
            subError = subError + (nodal(a) -pwLocal(a))^2 * w(a);
        end
    end
end

subError = sqrt(subError);
