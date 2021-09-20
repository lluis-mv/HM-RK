clc; clear all; close all;

% define some constants
syms h positive
syms M positive
syms k positive
syms dt positive
syms QBiot positive
syms AlphaStab real
syms DeltaPW positive

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
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + AlphaStab*h/12*[1, -1; -1 1];



for nNodes = [2:8, 15]
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
    
    % Estimate the stabilization factor
    for index = 2:nNodes-1
        Node1WaterPressure = NodalWaterPressure(index);
        Node2WaterPressure = NodalWaterPressure(index+1);
        disp([ 'Estimating the stabilizationFactor: mesh size ' int2str(nNodes), ' nodes. Between the nodes ' int2str(index) ])
        a = solve(Node1WaterPressure -Node2WaterPressure == 0, AlphaStab);
        disp(a)
    end
end

% Plot the obtained solution for several Gamma*AlphaStab
% First substitute the parameters

syms Gamma real;
hEval = 1/(nNodes-1); dtEval = 0.0001; MEval = 1; kEval=1; QBiotEval=1E9;
Solution = subs(Solution, h, hEval);
Solution = subs(Solution, dt, dtEval);
Solution = subs(Solution, M, MEval);
Solution = subs(Solution, k, kEval);
Solution = subs(Solution, QBiot, QBiotEval);
Solution = subs(Solution, DeltaPW, 1);

SavedSolution = Solution;

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
figure(1); clf;
plot(xx, pw, 'k', 'linewidth', 3); hold on
figure(2); clf;
for Gamma = 0:0.25:1.5
    AlphaEval = 3/MEval - 12 * dtEval * kEval /(hEval^2);
    if ( AlphaEval > 0)
        AlphaEval = Gamma * AlphaEval;
    else
        AlphaEval = 0;
    end
    Solution = subs(SavedSolution, AlphaEval);
    figure(1)
    plot( linspace(0,1,nNodes), 1-Solution(2:2:2*nNodes), '*-');
    xlabel('x/H'); ylabel('Water pressure')
    set(gca, 'FontSize', 12)
    hold on
    figure(2)
    plot( linspace(0,1,nNodes), Solution(1:2:2*nNodes-1), '*-');
    set(gca, 'FontSize', 12)
    xlabel('x/H'); ylabel('Displacement')
    hold on
end

figure(1)
legend('Analytical', '\gamma=0', '\gamma=0.25', '\gamma=0.5', '\gamma=0.75', '\gamma=1', '\gamma=1.25', '\gamma=1.5', 'location', 'best')
figure(2)
legend('\gamma=0', '\gamma=0.25', '\gamma=0.5', '\gamma=0.75', '\gamma=1', '\gamma=1.25', '\gamma=1.5', 'location', 'best')

figure(1)
print('../figures/MATLAB_FIGURES/water_pressure_1', '-dpdf')
figure(2)
print('../figures/MATLAB_FIGURES/displacement_1', '-dpdf')

