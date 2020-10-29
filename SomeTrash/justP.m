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



% Darcy law
H = -int( DN_DX*k*DN_DX', x, 0, h);

Mstab = AlphaStab*h/12*[1, -1; -1 1];

Mm = (1/M)*h/6*[2,1;1,2];


Ce = Mm+Mstab;
Ke = H;


DeltaPW = -1;

for nNodes = [2:8]%, 15]
    
    % Create the system matrix
    C = sym( zeros(nNodes, nNodes) );
    K = C;
    
    for i = 1:(nNodes-1)
        C( i-1 + [1:2], i-1 + [1:2]) = ...
            C( (i-1) + [1:2] , (i-1) + [1:2]) + Ce;
        K( (i-1) + [1:2], (i-1) + [1:2]) = ...
            K( (i-1) + [1:2] , (i-1) + [1:2]) + Ke;
    end
    
    % Apply dirichlet conditions
    % ZeroWaterPressure
    C(1, :) = 0;
    C(1,1) = 1;
    K(1,:) = 0;
    
    
    
    
    AMatrix = C\(K);
    
    Solution0 = sym(zeros(nNodes, 1));
    for i = 2:nNodes
        Solution0(i) = 1;
    end
    Solution0(1) = 0;
    
    
    A = C\(K);
    
    a = eig(A);
    aa = a(end);
    A2 = A*dt;
    a2 = eig(A2);
    aaa = a2(end)
    
    ii = eye(nNodes, nNodes);
    B = ii + dt * A;
    b = eig(B);
    bb = b(end);
    
    

    
    
    Solution = B*Solution0;
    
    NodalWaterPressure = Solution;
    
    
    % now Analytical solution
    Ba = subs(A, AlphaStab, 0);
    [vectors, values] = eig(Ba);
    x = 0*Solution;
    c = vectors\Solution0;
    for i = 1:size(values,1)
         x = x + c(i)*exp(values(i,i)*dt)*vectors(:,i);
    end
  
    difere = Solution-taylor(x, 'order', 3)
    
    
    % now second order
    B2 =  ii+0.5*dt*A + 0.5*dt*A*(ii+dt*A);
    Solution2 = B2*Solution0;
    difere2 = Solution2-taylor(x, 'order', 3)
    
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