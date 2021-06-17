function []=main2()

clc; clear all; close all;


i = 1;

NODES = [1, 2, 4, 8, 12, 16, 100, 200, 500,  1000];
nNodes = 20;
TEND = 0.001;
for nTimeSteps = NODES
    for StabMethod = [1:4]
        figure(3)
        subplot(4,1,StabMethod)
        [WP, U, L2, LInf] = ConstructAndSolveProblem(nNodes, TEND, nTimeSteps, StabMethod);
        L2Matrix(i,StabMethod) = L2;
        LInfMatrix(i,StabMethod) = LInf;
        hold on
        plot(linspace(0,1,nNodes), 1-WP)
        hold off
    end
    i = i+1
    pause(1)
end
figure(1)
%NODES = 1./NODES

NODES = TEND./NODES;
figure(1)       
for i = 1:4
    loglog(NODES, L2Matrix(:,i), '*-.')
    hold on
end
hold on
legend('LMV', 'LMV-1D', 'Borja', 'Primal')
ylabel('L2 Norm')


figure(2)       
for i = 1:4
    loglog(NODES, LInfMatrix(:,i), '*-.')
    hold on
end
hold on
legend('LMV', 'LMV-1D', 'Borja', 'Primal')
ylabel('Linf Norm')






function [NodalWaterPressure, NodalDisplacement, L2, LInf] = ConstructAndSolveProblem( nNodes, tFinal, nTimeSteps, StabMethod)

syms AlphaStab real

dt = tFinal/nTimeSteps;


nu =  0.0;
k =1;
DeltaPW = 1.0;

%QBiot = sym(1e12);

M =  1;
h = 1/(nNodes-1);
G = M*(1-2*nu)/2/(1-nu);
XX = linspace(0,1,nNodes)';

% Define shape functions (to then integrate)
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

% Define the elemental system
ElementalMatrix = sym(zeros(4,4));
ElementalMatrix2 = sym(zeros(4,4));

% Internal forces. Effective stress forces
ElementalMatrix([1, 3],[1, 3]) = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% Internal forces. Water pressure forces
ElementalMatrix([1, 3],[2, 4]) = -int( DN_DX * N', x, 0,  h);

% Mass conservation
% Volume change
ElementalMatrix([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
% Darcy law
ElementalMatrix([2,4],[2,4]) = dt*int( DN_DX*k*DN_DX', x, 0, h);
%Biot Coefficient
%ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + AlphaStab*h/12*[1, -1; -1 1];

% Matrix to multiply the old time
% Volume change
ElementalMatrix2([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
%Stabilization
ElementalMatrix2([2,4],[2,4]) = ElementalMatrix2([2,4],[2,4]) + AlphaStab*h/12*[1, -1; -1 1];


% Create the system matrix
SystemMatrix = sym( zeros(2*nNodes, 2*nNodes) );
SystemMatrix2 = sym( zeros(2*nNodes, 2*nNodes) );
    
for i = 1:(nNodes-1)
    SystemMatrix( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        SystemMatrix( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + ElementalMatrix;
    SystemMatrix2( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        SystemMatrix2( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + ElementalMatrix2;
end

% Apply dirichlet conditions
% ZeroWaterPressure
SystemMatrix(2, :) = 0;
SystemMatrix(2,2) = 1;

% Displacement
nn = 2*(nNodes -1)+1;
SystemMatrix(nn, :) = 0;
SystemMatrix(nn,nn) = 1;

% Apply dirichlet to the other matrix
SystemMatrix2(2,:) = 0;
SystemMatrix2(nn,:) = 0;


%Right hand side (water pressure at dirichlet)
FFExt = sym(zeros(2*nNodes, 1));
FFExt(2) = DeltaPW;

if (StabMethod == 1)
    AlphaEval = 3/M - 12 * dt * k / h^2;
    if (AlphaEval< 0)
        AlphaEval = 0;
    end
elseif ( StabMethod == 2)
    AlphaEval = 2/M - 12 * dt *k / h^2;
     if (AlphaEval< 0)
        AlphaEval = 0;
    end
elseif ( StabMethod == 3)
    AlphaEval = 1 / 2 / G;
else 
    AlphaEval = 0;
end
SystemMatrix = subs(SystemMatrix, AlphaStab, AlphaEval);
SystemMatrix = eval(SystemMatrix);

SystemMatrix2 = subs(SystemMatrix2, AlphaStab, AlphaEval);
SystemMatrix2 = eval(SystemMatrix2);

FFExt = eval(FFExt);


% Solve the linear systems
SolutionOld = zeros(2*nNodes, 1);
for t = 1:nTimeSteps
    Solution = SystemMatrix\(FFExt + SystemMatrix2 * SolutionOld);
    SolutionOld = Solution;
end


NodalWaterPressure = (zeros(nNodes, 1));
NodalDisplacement = (zeros(nNodes, 1));
for i = 1:nNodes
    NodalDisplacement(i) = Solution(2*(i-1)+1);
    NodalWaterPressure(i) = Solution(2*(i-1)+2);
end
    
        
% Compute Error Norms
TT =  M * tFinal*k;
nGP = 1;
[xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
xPlot = [];
yPlot = [];
for a = 1:size(xGP, 1)
    xPlot = [xPlot, xGP(a,:)];
    yPlot = [yPlot, pwGP(a,:)];
end
plot(xPlot, yPlot, 'r')
hold on
L2 = EvaluateL2( XX, pwGP, 1-NodalWaterPressure, nGP);

pw = EvalConsolidation(XX, TT);
plot(XX, pw, 'k')
LInf = max( abs(1-NodalWaterPressure - pw));




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
        nodal = mean(WPNodes(i:i+1));
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
