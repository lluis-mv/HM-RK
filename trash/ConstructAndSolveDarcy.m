%function [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)

function []=ConstructAndSolveDarcy()
i = 1;
close all; clf;
for nNodes = [10,20,50,60,70,80,100,200,500]
    [solution(i).t, solution(i).pw] = Aux(nNodes);
    pause(0.01)
    i = i+1;
end


function [timePointOfInterest, waterPressurePointOfInterest] = Aux(nNodes)

H = 1;

%nNodes = 1000;
tFinal = 1e-4;
nTimeSteps = 5000;
nPlot = 1000;

StabMethod = 1;
nGP = 3;

XX = linspace(0,H,nNodes)';

xPointOfInterest = 0.65*H;
nPointOfInterest = find( min(abs(XX-xPointOfInterest))==abs(XX-xPointOfInterest));
nPointOfInterest = nPointOfInterest(1);
timePointOfInterest = zeros(1,nTimeSteps);
timePointOfInterest = [1:nTimeSteps]/nTimeSteps*tFinal;
waterPressurePointOfInterest = zeros(1,nTimeSteps);
dt = tFinal/nTimeSteps;
QBiot = 100;
OneOverBiotQ = 1/QBiot;

beta = 0.3025;
gamma = 0.6;

alpha1 = 1 / beta / (dt^2);
alpha2 = 1 / beta / dt;
alpha3 = 1/2/beta -1;
alpha4 = gamma / beta/dt;
alpha5 = 1 - gamma/beta;
alpha6 = (1-gamma/beta/2)*dt;

E = 500e4;
nu = 0.0;
M =  E*(1-nu)/(1+nu)/(1-2*nu);
QBiot = M;

OneOverBiotQ = 0.00001/QBiot;
OneOverBiotQ = 1/QBiot;
OneOverBiotQ = 1e-3;
k = 1e-3;


rho_mixture =0;
rho_water = 0;
porosity = 0.4;
rho_mixture = rho_mixture*(1-porosity) + porosity * rho_water;

% if ( nargin < 5)
%     Draw = false;
% end
% if ( nargin < 6)
%     nu = 1/3;
% end


h = H/(nNodes-1);
G = M*(1-2*nu)/2/(1-nu);


AlphaStab = 0;
%AlphaStab = 2/M+2*OneOverBiotQ - beta/gamma * 12 * dt *k / h^2;
%AlphaStab = 2/M + 2*OneOverBiotQ;
%AlphaStab = AlphaStab*(0.5 + 0.5*tanh( 2 - 12*dt * k * M/h^2) );
if (AlphaStab < 0)
    AlphaStab = 0;
end



% Define shape functions (to then integrate)4
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

StiffnessMatrix = zeros(6,6);
DampingMatrix = zeros(6,6);
MassMatrix = zeros(6,6);
% Internal forces
StiffnessMatrix([1, 4],[1, 4]) = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% water forces
StiffnessMatrix([1, 4],[3, 6]) = -int( DN_DX * N', x, 0,  h);


%Darcy law
StiffnessMatrix([2, 5],[3 ,6]) = +int( N*DN_DX', x, 0, h);
DampingMatrix([2, 5],[2 ,5]) = +int( (1/k)*N*N', x, 0, h);

%Mass Conservation
DampingMatrix([3,6],[1,4]) = -int( N*DN_DX', x, 0, h);
DampingMatrix([3,6],[2,5]) = -int( N*DN_DX', x, 0, h);


%Biot Coefficient
DampingMatrix([3,6],[3,6]) = -(OneOverBiotQ)*int( N*N', x, 0, h);
%Stabilization
DampingMatrix([3,6],[3,6]) = DampingMatrix([3,6],[3,6]) - AlphaStab*h/12*[1, -1; -1 1];

% dynamic terms;
MassMatrix([1,4],[1,4]) = rho_mixture * int( N*N', x, 0, h);
MassMatrix([2,5],[1,4]) = rho_water * int( N*N', x, 0, h);
MassMatrix([1,4],[2,5]) = rho_water * int( N*N', x, 0, h);
MassMatrix([2,5],[2,5]) = rho_water/porosity * int( N*N', x, 0, h);


% Create the system matrix
SystemMatrix = sparse( 3*nNodes, 3*nNodes );
DampingMatrixBig = SystemMatrix;
MassMatrixBig = DampingMatrixBig;
Solution= zeros(3*nNodes,1);
Solution([3:3:3*nNodes]) = 0;
SolutionOld = Solution;
Velocity  = zeros(3*nNodes,1);
VelocityOld = zeros(3*nNodes,1);
Acceleration = Velocity;
AccelerationOld = Velocity;

StiffnessMatrix = (StiffnessMatrix);
DampingMatrix =(DampingMatrix);
MassMatrix = (MassMatrix);

KMatrix = alpha1 * MassMatrix + alpha4 * DampingMatrix + StiffnessMatrix;

% Before solving the problem I should  put Now and
for i = 1:(nNodes-1)
    SystemMatrix( 3*(i-1) + [1:6],3*(i-1) + [1:6]) = ...
        SystemMatrix( 3*(i-1) + [1:6],3*(i-1) + [1:6]) + KMatrix;
    DampingMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) = ...
        DampingMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) + DampingMatrix;
    MassMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) = ...
        MassMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) + MassMatrix;
end

% Apply dirichlet conditions
% ZeroWaterPressure
% nnWP = [3,3*nNodes];
% BC1
nnWP = [3];
nnDW = 3*(nNodes -1)+2;
%BC2
nnWP =[3,3*nNodes];
nnDW = [];

SystemMatrix(nnWP, :) = 0;
SystemMatrix(nnWP,nnWP) = eye(length(nnWP));

% Displacement
nnD = 1:3:3*nNodes;
SystemMatrix(nnD, :) = 0;
SystemMatrix(nnD,nnD) = eye(length(nnD));

% Displacement water
SystemMatrix(nnDW, :) = 0;
SystemMatrix(nnDW,nnDW) = eye(length(nnDW));




NodalWaterPressure = zeros(nNodes, 1);
NodalDisplacement = NodalWaterPressure;
NodalWaterDisplacement = NodalWaterPressure;

FFExt = sparse(3*nNodes, 1);
FFExt(3) = 0;
FFExt(1) = 0;
thisNode = find( min( abs(XX-xPointOfInterest)) == abs(XX-xPointOfInterest))
thisNode = thisNode(1);
thisNode = [thisNode, thisNode+1];
% FFExt(3*(thisNode-1)+3) = 1.00/h;


FF = FFExt;
Solution(nnWP) = zeros(length(nnWP),1);
SolutionOld(nnWP) = zeros(length(nnWP),1);
FFExt(nnWP) = 0;
FF(nnWP) = 0;
FFExt(end) = 0.1;
SolutionOld(end) = FFExt(end);
for step = 1:nTimeSteps
    FF = FFExt + DampingMatrixBig * ( alpha4 * SolutionOld - alpha5*VelocityOld - alpha6 * AccelerationOld);
    FF = FF +  MassMatrixBig * ( alpha1 * SolutionOld + alpha2*VelocityOld + alpha3 * AccelerationOld);
    FF(nnWP) = FFExt(nnWP);
    FF(nnD) = FFExt(nnD);
    FF(nnDW) = FFExt(nnDW);
    Solution = SystemMatrix\FF;
    dSolution = Solution-SolutionOld;
    Acceleration = alpha1 * dSolution - alpha2 * VelocityOld - alpha3 * AccelerationOld;
    Velocity = alpha4*dSolution + alpha5 * VelocityOld + alpha6 * AccelerationOld;
    AccelerationOld = Acceleration;
    VelocityOld = Velocity;
    SolutionOld = Solution;
    
    for i = 1:nNodes
        NodalDisplacement(i) = Solution(3*(i-1)+1);
        NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
        NodalWaterPressure(i) = Solution(3*(i-1)+3);
    end
    timePointOfInterest(step) = step/nTimeSteps*tFinal;
    waterPressurePointOfInterest(step) = -NodalWaterPressure(nPointOfInterest);
    if (step/nPlot == floor(step/nPlot) || step == 1)
        TT =  M *k*step/nTimeSteps * tFinal;
        figure(1);
        subplot(3,2,1)
        %uu = EvalConsolidationDISPL(XX, TT);
        %plot(XX, uu, 'k')
        %hold on
        plot(XX, NodalDisplacement)
        hold off
        subplot(2,2,2)
        plot(XX, NodalWaterDisplacement)
        subplot(2,2,3)
        plot(XX, -NodalWaterPressure)
        
        subplot(2,2,4)
        plot(XX, Velocity(2:3:3*nNodes))
        
        
        pause(0.001)
    end
end






