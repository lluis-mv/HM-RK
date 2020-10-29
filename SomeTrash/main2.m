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
    % Solve the system of matrices
    Solution = Solution0 + dt* ( -C\(K*Solution0) )
    
    %Solution = Solution + dt* ( -C\(K*Solution))
    %Solution = Solution + dt* ( C\(K*Solution))
    
    
    
    % Other way, LMV
    aa = eye(2*nNodes, 2*nNodes) - (dt*AMatrix);
    bb = eig(aa);
    m = bb(end);
    
    
    m = subs(m, h, 0.01)
    m = subs(m, M, 1);
    m = subs(m, k, 1);
    m = subs(m, AlphaStab, 0);
    ezplot(m);
    
    m2 = bb(end);
    m2 = subs(m2, h, 0.01)
    m2 = subs(m2, M, 1);
    m2 = subs(m2, k, 1);
    AA = solve(m2+1, AlphaStab)
    m2 = subs(m2, AlphaStab, AA);
    
    hold on
    ezplot(m2);
    hold off
    
    
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
    
    
    % Now the same but second order
    
    K1 = ( -C\(K*Solution0) );
    
    X2 = Solution0 + dt*1*K1;
    K2 = ( -C\(K*X2) );
    
    XEnd = Solution0 + 0.5*dt*K1 + 0.5*K2*dt;
    
    
    ii = eye(2*nNodes, 2*nNodes)
    A = -(C\K); 
    B = (ii + 0.5*dt*A + 0.5*dt*A*(ii+dt*A));
    XEnd2 = B*Solution0;
    
    
    
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