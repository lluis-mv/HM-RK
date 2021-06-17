%function [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)

function []=ConstructAndSolveProblemUWwP()
i = 1;
close all; clf;
for nNodes = [20,50,100,200,500]
    [solution(i).t, solution(i).pw] = Aux(nNodes);
    plot(solution(i).t, solution(i).pw);
    hold on
    pause(0.01)
    i = i+1;
end
legend('20','50','100','200','500','location','best')

function [timePointOfInterest, waterPressurePointOfInterest] = Aux(nNodes)

H = 1;

%nNodes = 1000;
tFinal = 5e-4;
nTimeSteps = 5000;
StabMethod = 1;
nPlot = 100000;
nGP = 3;

XX = linspace(0,H,nNodes)';

xPointOfInterest = 0.2;
nPointOfInterest = find( min(abs(XX-xPointOfInterest))==abs(XX-xPointOfInterest));
nPointOfInterest = nPointOfInterest(1);
timePointOfInterest = zeros(1,nTimeSteps);
timePointOfInterest = [1:nTimeSteps]/nTimeSteps*tFinal;
waterPressurePointOfInterest = zeros(1,nTimeSteps);
dt = tFinal/nTimeSteps;
QBiot = 1;
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
OneOverBiotQ = 1/QBiot;

k = 0.001/10;
rho_mixture =2.650;
rho_water = 1.0;
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

DN_DX = eval(DN_DX)
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
nnWP = 3;
SystemMatrix(nnWP, :) = 0;
SystemMatrix(nnWP,nnWP) = 1;

% Displacement
nnD = 3*(nNodes -1)+1;
SystemMatrix(nnD, :) = 0;
SystemMatrix(nnD,nnD) = 1;

% Displacement water
nnDW = 3*(nNodes -1)+2;

SystemMatrix(nnDW, :) = 0;
SystemMatrix(nnDW,nnDW) = 1;

%Inverse = inv(SystemMatrix);

NodalWaterPressure = zeros(nNodes, 1);
NodalDisplacement = NodalWaterPressure;
NodalWaterDisplacement = NodalWaterPressure;

FFExt = sparse(3*nNodes, 1);
FFExt(3) = -1;
FFExt(1) = -1;
FF = FFExt;
Solution(nnWP) = FFExt(3);
SolutionOld(nnWP) = FFExt(3);
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
    if (step/nPlot == floor(step/nPlot))
        TT =  M *k*step/nTimeSteps * tFinal;
        figure(1);
        subplot(3,2,1)
        %uu = EvalConsolidationDISPL(XX, TT);
        %plot(XX, uu, 'k')
        %hold on
        plot(XX, NodalDisplacement)
        hold off
        subplot(3,2,2)
        plot(XX, NodalWaterDisplacement)
        subplot(3,2,3)
        %pw = EvalConsolidation(XX, TT);
        %plot(XX, pw, 'k')
        %hold on
        plot(XX, -NodalWaterPressure)
        ylim([0,2])
        xlim([0,H])
        hold off
        pause(0.01)
        %TryToComputeSigma
        xx = 2*nNodes-1;
        yy = 2*nNodes-1;
        for i = 1:nNodes-1
            index = [2*i-1,2*i];
            xx(index) = XX([i,i+1]);
            yy(index) = M*DN_DX'*NodalDisplacement([i,i+1]);
        end
        subplot(3,2,4)
        plot(xx,yy)
        ylim([0,2])
        xlim([0,H])
        
        xx = 2*nNodes-1;
        yy = 2*nNodes-1;
        for i = 1:nNodes-1
            index = [2*i-1,2*i];
            xx(index) = XX([i,i+1]);
            yy(index) = M*DN_DX'*NodalDisplacement([i,i+1])-mean(NodalWaterPressure([i,i+1]));
        end
        subplot(3,2,5)
        plot(xx,yy)
        ylim([0,2])
        xlim([0,H])
        subplot(3,2,6)
        figure(3)
        subplot(1,1,1)
        plot(timePointOfInterest, waterPressurePointOfInterest)
        grid minor
        pause(0.01)
        figure(1)
        pause(0.01)
    end
end






