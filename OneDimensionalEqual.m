function []=OneDimensionalEqual()

nodes = linspace(0, 1, 200);


nSteps = 200;
RKMethod = 1;

M = 1;
k = 1;


substepping = false;
if ( RKMethod > 10)
    alfa = floor(RKMethod/10);
    substepping = true;
end



dt = 1E-1;
dt = dt/nSteps;
h = nodes(2)-nodes(1);
AlphaStab = 12*M*dt*k/(M*h^2)
AlphaStab = dt*(-3*h^2 + 10*M*k)/(M*h^2)
if (substepping)
    dt2 = dt/alfa;
    AlphaStab = dt2*(-3*h^2 + 8*M*k)/(M*h^2)
end
% AlphaStab = 0;
% AlphaStab = -(3*h^2 + 12*M*k)/(M*h^2)

% AlphaStab = ((12*M*dt*k)/h^2 - 3)/M
if (AlphaStab < 0)
    AlphaStab = 0;
end
AlphaStab


[AMatrix, FFExt, X0] = ConstructMatrices(nodes, M, k, AlphaStab, dt);

[X, NodalWaterPressure] = IntegrateProblem(AMatrix, FFExt, X0, RKMethod, dt, nSteps);


[x, y] = EvaluateSolution(nSteps*dt, M*k);
plot(x, y, 'k-.', 'linewidth', 2.5); hold on
plot(nodes, NodalWaterPressure, '*-')
hold off




function [X, NodalWP] = IntegrateProblem(AMatrix, FFExt, X, RKMethod, dt, nSteps)

[aRK, bRK] = GetRungeKutta(RKMethod);

k = zeros( length(X), length(bRK));

for t = 1:nSteps
    for i = 1:length(bRK)
        Xstep = X;
        for j = 1:i-1
            Xstep = Xstep + dt*aRK(i,j)*k(:,j);
        end
        
        [k(:,i)] = -AMatrix*Xstep+FFExt;
        
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
    for i = 1:nNodes
        NodalWP(i) = X(2*(i-1)+2);
    end

% % 
% % 
% % for tt = 1:nSteps
% %     Solution = Solution0 + dt* (-AMatrix*Solution0 + FFExt );
% %     if (tt > 0)
% %         FFExt = 0*FFExt;
% %     end
% %     
% %     Solution0 = Solution;
% %     
% %     NodalWaterPressure = (zeros(nNodes, 1));
% %     NodalDisplacement = (zeros(nNodes, 1));
% %     for i = 1:nNodes
% %         NodalDisplacement(i) = Solution(2*(i-1)+1);
% %         NodalWaterPressure(i) = Solution(2*(i-1)+2);
% %     end
% %     
% % %     plot(X, NodalWaterPressure, '*-'); hold on
% % %     pause(0.01)
% % end




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
AMatrix(2,:) = 0;

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
Ke = Ke(ind,ind);