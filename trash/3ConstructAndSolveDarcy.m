%function [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)

function []=v3ConstructAndSolveDarcy()
i = 1;
close all; clf;
for nNodes = [10,20,50,60,70,80,100,200,500]
    [solution(i).t, solution(i).pw] = Aux(nNodes);
    pause(0.01)
    i = i+1;
end


function [timePointOfInterest, waterPressurePointOfInterest] = Aux(nNodes)

OneOverBiotQ = 1/QBiot;
k = 1e-3;

H = 1;

%nNodes = 1000;
tFinal = 1e-4;
nTimeSteps = 5000;
nPlot = 1000;

StabMethod = 1;


XX = linspace(0,H,nNodes)';

dt = tFinal/nTimeSteps;
QBiot = 100;





h = H/(nNodes-1);


AlphaStab = 0;
if (AlphaStab < 0)
    AlphaStab = 0;
end



% Define shape functions (to then integrate)4
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

StiffnessMatrix = zeros(6,6);
DampingMatrix = zeros(6,6);

% Internal forces
StiffnessMatrix([1, 4],[1, 4]) = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% water forces
StiffnessMatrix([1, 4],[3, 6]) = -int( DN_DX * N', x, 0,  h);


%Darcy law
StiffnessMatrix([2, 5],[3 ,6]) = +int( N*DN_DX', x, 0, h);
StiffnessMatrix([2, 5],[2 ,5]) = (1/k)*int( N*N', x, 0, h);

StiffnessMatrix([3,6],[2,5]) = +int( N*DN_DX', x, 0, h);


%Biot Coefficient
DampingMatrix([3,6],[3,6]) = +(OneOverBiotQ)*int( N*N', x, 0, h);



% Create the system matrix
SystemMatrix = sparse( 3*nNodes, 3*nNodes );
DampingMatrixBig = SystemMatrix;

Solution= zeros(3*nNodes,1);
Solution([3:3:3*nNodes]) = 0;
SolutionOld = Solution;



% Before solving the problem I should  put Now and
for i = 1:(nNodes-1)
    SystemMatrix( 3*(i-1) + [1:6],3*(i-1) + [1:6]) = ...
        SystemMatrix( 3*(i-1) + [1:6],3*(i-1) + [1:6]) + StiffnessMatrix + 1/dt*DampingMatrix;
    DampingMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) = ...
        DampingMatrixBig( 3*(i-1) + [1:6],3*(i-1) + [1:6]) + 1/dt*DampingMatrix;
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

SolutionOld(end) = FFExt(end);
for step = 1:nTimeSteps
    Solution = SystemMatrix\(FF+DampingMatrix*SolutionOld);
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






