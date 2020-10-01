%clc; 
clear all; close all;

% define some constants
syms h positive
syms M positive
syms k positive
syms dt positive
syms QBiot positive
syms AlphaStab real
% syms DeltaPW positive

% Define shape functions (to then integrate)
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

% Define the elemental system


% Internal forces. Effective stress forces
K = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% Internal forces. Water pressure forces
Qt = -int( DN_DX * N', x, 0,  h);


% Mass conservation
% Volume change
Q = int( N*DN_DX', x, 0, h);
% Darcy law
H = int( DN_DX*k*DN_DX', x, 0, h);

Mstab = AlphaStab*h/12*[1, -1; -1 1];




Ce = [K, Qt; Q, Mstab];
Ke = [zeros(2,4); zeros(2,2), H];

ind = [1,3,2,4];
Ce = Ce(ind,ind);
Ke = Ke(ind,ind);


DeltaPW = -1;

for nNodes = [2:8]%, 15]
    % Create the system matrix
    C = sym( zeros(2*nNodes, 2*nNodes) );
    K = C;
    
    for i = 1:(nNodes-1)
        C( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
            C( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + Ce;
        K( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
            K( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + Ke;
    end
    
    % Apply dirichlet conditions
    % ZeroWaterPressure
    C(2, :) = 0;
    C(2,2) = 1;
    K(2,:) = 0;
    
    % Displacement
    nn = 2*(nNodes -1)+1;
    C(nn, :) = 0;
    C(nn,nn) = 1;
    
    %Right hand side (water pressure at dirichlet)
    FFExt = sym(zeros(2*nNodes, 1));
    FFExt(2) = DeltaPW/dt;
    
    Solution0 = 0*FFExt;
    for i = 1:nNodes
        Solution0(2*i) = 1;
    end
    Solution0(2) = 0;
    % Solve the system of matrices
    Solution = Solution0 + dt* ( -C\(K*Solution0) + 0*FFExt)
    
    Solution = Solution + dt* ( -C\(K*Solution))
    %Solution = Solution + dt* ( C\(K*Solution))
    
    
    % Other way, LMV
    aa = eye(2*nNodes, 2*nNodes) - dt*C\K;
    bb = eig(aa);
    
    if (nNodes > 2)
        e1 = bb(end);
        e2 = bb(end-1);
        
        a1 = solve(e1+1, AlphaStab);
        a1 = a1(1);
        a2 = solve(e2+1, AlphaStab);
        a2 = a2(1);
        
        MM = solve(a2, M);
        a23 = subs(a1, M, MM);
        simplify(a23)
    end
    
    for i = 1:size(aa,2)
        term = aa(i, end)
        aaa = solve(term, AlphaStab)
        disp(aaa)
        hola = 1;
    end
    
    NodalWaterPressure = sym(zeros(nNodes, 1));
    NodalDisplacement = sym(zeros(nNodes, 1));
    for i = 1:nNodes
        NodalDisplacement(i) = Solution(2*(i-1)+1);
        NodalWaterPressure(i) = Solution(2*(i-1)+2);
    end
  
    hola = 1;
    NodalWaterPressure = simplify(NodalWaterPressure);
    
    % Estimate the stabilization factor
    for index = 1:nNodes-1
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
