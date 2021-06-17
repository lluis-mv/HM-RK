%function [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)

function []=v3ConstructAndSolveDarcy()
i = 1;
close all; clf;
for nNodes = [10,20,50,60,70,80,100,200,500]
    Aux(nNodes);
    pause(0.01)
    i = i+1;
end


function [MySystemMatrix] = Aux(nNodes)

OneOverBiotQ = 1e-6;
% OneOverBiotQ = 0;
k = 1e-6;

H = 0.1;

%nNodes = 1000;
tFinal = 1e-8;
nTimeSteps = 100;
nPlot = 100;



XX = linspace(0,H,nNodes)';

dt = tFinal/nTimeSteps;

h = H/(nNodes-1);

denom = (k)
nume = OneOverBiotQ/dt;

ratio = denom/nume





AlphaStab = 5*OneOverBiotQ/dt;
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
StiffnessMatrix([1, 4],[1, 4]) = +int(  DN_DX * (DN_DX'), x, 0,  h);
% water forces
StiffnessMatrix([1, 4],[3, 6]) = -int( DN_DX * N', x, 0,  h);


%Darcy law
StiffnessMatrix([2, 5],[3 ,6]) = +int( N*DN_DX', x, 0, h);
StiffnessMatrix([2, 5],[2 ,5]) = (1/k)*int( N*N', x, 0, h);

StiffnessMatrix([3,6],[2,5]) = +int( N*DN_DX', x, 0, h);


%Biot Coefficient
DampingMatrix([3,6],[3,6]) = +(OneOverBiotQ)*int( N*N', x, 0, h);

DampingMatrix([3,6],[3,6]) = DampingMatrix([3,6],[3,6]) + AlphaStab*h/12*[1, -1; -1 1];

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
%BC3
nnWP = [];
nnDW = [2,3*nNodes-1]

SystemMatrix(nnWP, :) = 0;
SystemMatrix(nnWP,nnWP) = eye(length(nnWP));

% Displacement
nnD = 1:3:3*nNodes;
SystemMatrix(nnD, :) = 0;
SystemMatrix(nnD,nnD) = eye(length(nnD));

% Displacement water
SystemMatrix(nnDW, :) = 0;
SystemMatrix(nnDW,nnDW) = eye(length(nnDW));

DampingMatrixBig(nnWP,:) = 0;


NodalWaterPressure = zeros(nNodes, 1);
NodalDisplacement = NodalWaterPressure;
NodalWaterDisplacement = NodalWaterPressure;

FFExt = sparse(3*nNodes, 1);



FF = FFExt;
Solution(nnWP) = zeros(length(nnWP),1);
SolutionOld(nnWP) = zeros(length(nnWP),1);
FFExt(nnWP) = 0;
FF(nnWP) = 0;
FF(end) = 0;
SolutionOld(end) = FFExt(end);
% FF(nNodes*3/2) = 0.01;
% FF(nNodes*3/2+3) = 0.02;
% FF(nNodes*3/2+6) = 0.012;
% FF(nNodes*3/2+9) = 0.005;
FF(2) = 0.01;
for step = 1:nTimeSteps
    Solution = SystemMatrix\(FF+DampingMatrixBig*SolutionOld);
    SolutionOld = Solution;
    indexes = [2:3:3*nNodes, 3:3:3*nNodes];
    MySystemMatrix = SystemMatrix(indexes, indexes);
    for i = 1:nNodes
        NodalDisplacement(i) = Solution(3*(i-1)+1);
        NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
        NodalWaterPressure(i) = Solution(3*(i-1)+3);
    end
    
    if (step/nPlot == floor(step/nPlot) || step == 1)
    
        figure(1);
        subplot(3,2,1)
        plot(XX, NodalWaterDisplacement)

        subplot(2,2,2)
        plot(XX, NodalWaterDisplacement)
        subplot(2,2,3)
        plot(XX, -NodalWaterPressure)
        

        pause(1)
    end
end






