function []=main()

clc; clear all; close all;

% define some constants
syms h positive
syms M positive
syms k positive
syms dt positive
syms QBiot positive
syms AlphaStab real
syms DeltaPW positive

nGP = 2;

nNodesD = 5;
E = 100.0;
nu =  0.3;
k = sym(1e-6);
DeltaPW = 1.0;
QBiot = sym(1e12);
M =  E*(1-nu)/ (1+nu) / (1-2*nu);
h = 1/(nNodesD-1);
G = E/2/(1+nu);
XX = linspace(0,1,nNodesD)';

% Define shape functions (to then integrate)
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

% Define the elemental system

ElementalMatrix = sym(zeros(4,4));
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



for nNodes = [nNodesD]
    % Create the system matrix
    SystemMatrix = sym( zeros(2*nNodes, 2*nNodes) );
    
    for i = 1:(nNodes-1)
        SystemMatrix( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
            SystemMatrix( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + ElementalMatrix;
    end
    
    % Apply dirichlet conditions
    % ZeroWaterPressure
    SystemMatrix(2, :) = 0;
    SystemMatrix(2,2) = 1;
    
    % Displacement
    nn = 2*(nNodes -1)+1;
    SystemMatrix(nn, :) = 0;
    SystemMatrix(nn,nn) = 1;
    
    %Right hand side (water pressure at dirichlet)
    FFExt = sym(zeros(2*nNodes, 1));
    FFExt(2) = DeltaPW;
    
    % Solve the system of matrices
    Solution = SystemMatrix\FFExt;
    
    
    NodalWaterPressure = sym(zeros(nNodes, 1));
    NodalDisplacement = sym(zeros(nNodes, 1));
    for i = 1:nNodes
        NodalDisplacement(i) = Solution(2*(i-1)+1);
        NodalWaterPressure(i) = Solution(2*(i-1)+2);
    end
    
    NodalWaterPressure = simplify(NodalWaterPressure);
    
end

% Plot the obtained solution for several Gamma*AlphaStab
% First substitute the parameters


SavedSolution = Solution;

% plot analytical solution
abc = linspace(-3,6, 20);

ddTT = 10.^abc;
ErrorInfiPrimal = [];
ErrorInfiJo = [];
ErrorInfiBorja = [];

ErrorL2Primal = [];
ErrorL2Jo = [];
ErrorL2Borja = [];

index = 1;

for dtEval = ddTT
    disp(dtEval)
    pause(0.0001)
    
    TT = M * dtEval*k;
    pw = EvalConsolidation(XX, TT);
    
    [xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
    
    WPPrimal = 1-NodalWaterPressure;
    WPPrimalDT = subs(WPPrimal, dt, dtEval);
    WPPrimal = subs(WPPrimalDT, AlphaStab, 0);
    WPPrimal = eval(WPPrimal);
    subError = max( abs(WPPrimal - pw));
    ErrorInfiPrimal = [ErrorInfiPrimal, subError];
    subError = EvaluateL2( XX, pwGP, WPPrimal, nGP);
    ErrorL2Primal = [ErrorL2Primal, subError];
    
    
    AlphaEval = 3/M - 12 * dtEval * k / h^2;
    if (AlphaEval< 0)
        AlphaEval = 0;
    end
    WPPrimal = subs(WPPrimalDT, AlphaStab, AlphaEval);
    WPPrimal = eval(WPPrimal);
    subError = max( abs(WPPrimal - pw));
    ErrorInfiJo = [ErrorInfiJo, subError ];
    subError = EvaluateL2( XX, pwGP, WPPrimal, nGP);
    ErrorL2Jo = [ErrorL2Jo, subError];
    
    AlphaEval = 1/2/G;
    WPPrimal = subs(WPPrimalDT, AlphaStab, AlphaEval);
    WPPrimal = eval(WPPrimal);
    subError = max( abs(WPPrimal - pw));
    ErrorInfiBorja = [ErrorInfiBorja, subError ];
    subError = EvaluateL2( XX, pwGP, WPPrimal, nGP);
    ErrorL2Borja = [ErrorL2Borja, subError];
    
    figure(1)
    loglog( ddTT(1:index), ErrorInfiPrimal, 'ks:')
    hold on
    loglog( ddTT(1:index), ErrorInfiJo, 'r>-.')
    loglog( ddTT(1:index), ErrorInfiBorja, 'g*:')
    crit = 3 / M * h^2 / k;
    loglog(crit*[1,1], [1e-2, 1], 'c')
    hold off
    pause(0.01)
    
    figure(2)
    loglog( ddTT(1:index), ErrorL2Primal, 'ks:')
    hold on
    loglog( ddTT(1:index), ErrorL2Jo, 'r>-.')
    loglog( ddTT(1:index), ErrorL2Borja, 'g*:')
    pause(0.01)
    crit = 3 / M * h^2 / k;
    loglog(crit*[1,1], [1e-2, 1], 'c')
    hold off
    pause(0.01)
    
    index = index +1
end





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
    N = [0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3);
        0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3)];
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
    N = [0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3);
        0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3)];
    L = XX(2)-XX(1);
    w = L/2;
    for i = 1:length(XX)-1
        nodal = N * [WPNodes(i:i+1)];
        pwLocal = pwGP(i,:);
        for a = 1:2
            subError = subError + (nodal(a) -pwLocal(a))^2 * w;
        end
    end
end

subError = sqrt(subError);
