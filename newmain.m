%clc;
clear all; close all;

% define some constants
syms h real
syms M real
syms k real
syms dt real
syms QBiot real
syms AlphaStab real


% Define shape functions (to then integrate)
syms x real
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
    
    
    AMatrix = C\(-K);
    
    Solution0 = sym(zeros(2*nNodes, 1));
    for i = 1:nNodes
        Solution0(2*i) = 1;
    end
    Solution0(2) = 0;
    
    
    A = C\(-K);
    A = simplify(A);
    
    a = eig(A);
    aa = a(end);
    A2 = A*dt;
    a2 = eig(A2);
    aaa = a2(end)
    
    B = eye(2*nNodes, 2*nNodes) + dt * A;
    b = eig(B);
    bb = b(end);
    
    ii = eye(2*nNodes, 2*nNodes);
    B2 =  ii+0.5*dt*A + 0.5*dt*A*(ii+dt*A);
%     B2 = ii + 0.5* dt * A + 0.5*A*(ii+dt*A);
    b2 = eig(B2)
    bb2 = b2(end);
    
    Solution = B*Solution0;
    
    
    % now Analytical solution
    Ba = subs(A, AlphaStab, 0);
    [vectors, values] = eig(Ba);
    x = 0*Solution;
    c = vectors\Solution0;
    for i = 1:size(values,1)
        x = x + c(i)*exp(values(i,i)*dt)*vectors(:,i);
    end
    
    difere = Solution-taylor(x, 'order', 3)
    
    difere = subs(difere, AlphaStab, 12*dt*k/h^2 + 3/M)
    
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