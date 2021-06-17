%function [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)
clear all;
% OneOverBiotQ = 1e-6;
% % OneOverBiotQ = 0;
% k = 1e-6;
% 
% H = 0.1;
% 
% %nNodes = 1000;
% tFinal = 1e-8;
% nTimeSteps = 100;
% nPlot = 100;
% 
% 
% 
% XX = linspace(0,H,nNodes)';
% 
% dt = tFinal/nTimeSteps;
% 
% h = H/(nNodes-1);
% 
% denom = (k)
% nume = OneOverBiotQ/dt;
% 
% ratio = denom/nume


syms OneOverBiotQ positive;
syms k positive;
syms h positive;
syms dt positive;
% dt = 1e-7
syms AlphaStab real;

nNodes =4;

% OneOverBiotQ = 0.00001;

% Define shape functions (to then integrate)4
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

StiffnessMatrix = sym(zeros(6,6));
DampingMatrix = sym(zeros(6,6));

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
SystemMatrix = sym(zeros( 3*nNodes, 3*nNodes ));
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
nnDW = [2,3*nNodes-1];
% nnDW = 2:3:3*nNodes;
nnDW = [2];
nnWP = 3*nNodes;
if ( length(nnWP) > 0)
    SystemMatrix(nnWP, :) = 0;
    SystemMatrix(nnWP,nnWP) = eye(length(nnWP));
end

% Displacement
nnD = 1:3:3*nNodes;
SystemMatrix(nnD, :) = 0;
SystemMatrix(nnD,nnD) = eye(length(nnD));

% Displacement water
SystemMatrix(nnDW, :) = 0;
SystemMatrix(nnDW,nnDW) = eye(length(nnDW));

DampingMatrixBig(nnWP,:) = 0;


% NodalWaterPressure = zeros(nNodes, 1);
% NodalDisplacement = NodalWaterPressure;
% NodalWaterDisplacement = NodalWaterPressure;

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
for step = 1:1
    Solution = SystemMatrix\(FF+DampingMatrixBig*SolutionOld);
    SolutionOld = Solution;
    indexes = [2:3:3*nNodes, 3:3:3*nNodes];
    MySystemMatrix = SystemMatrix(indexes, indexes);
    
    for i = 1:nNodes
        NodalDisplacement(i) = Solution(3*(i-1)+1);
        NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
        NodalWaterPressure(i) = Solution(3*(i-1)+3);
    end
  
    n1 = NodalWaterDisplacement(2);
    n2 = NodalWaterDisplacement(1);
    difere = n1-n2
    al1 = solve(difere, AlphaStab)
    
    n1 = NodalWaterPressure(1);
    n2 = NodalWaterPressure(2);
    difere = n1-n2
    al2 = solve(difere, AlphaStab)
end






