
function [ProblemData] = ConstructAndSolveProblem(ProblemData)


nNodes = ProblemData.nNodes;
tFinal = ProblemData.tFinal;
nTimeSteps = ProblemData.nTimeSteps;
StabMethod = ProblemData.StabMethod;
k = ProblemData.k;
M = ProblemData.M;
nu = ProblemData.nu;
QBiot = ProblemData.QBiot;


dt = tFinal/nTimeSteps;

h = 1/(nNodes-1);
XX = linspace(0,1,nNodes)';

if (StabMethod == 1)
    AlphaEval = 0;
elseif ( StabMethod == 2)
    G = M*(1-2*nu)/2/(1-nu);
    AlphaEval = 1 / 2 / G;
elseif (StabMethod == 3)
    AlphaEval = 2/M - 6 * dt * k / h^2;
    AlphaEval = AlphaEval * (0.5 + 0.5*tanh( 2 - 12*dt * k * M/h^2) );
elseif ( StabMethod == 4)
    AlphaEval = 3/M;
elseif ( StabMethod == 5)
    AlphaEval = 2/M - 12 * dt * k / h^2;
    if (AlphaEval< 0)
        AlphaEval = 0;
    end
else
    AlphaEval = 3/M - 12 * dt * k / h^2;
    if (AlphaEval< 0)
        AlphaEval = 0;
    end
end

if ( AlphaEval < -1e-8)
    AlphaEval = 0;
end

% Define shape functions (to then integrate)
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

% Define the elemental system
ElementalMatrix = sym(zeros(4,4));
ElementalMatrix2 = sym(zeros(4,4));

% Internal forces. Effective stress forces
ElementalMatrix([1, 3],[1, 3]) = int( M *  DN_DX * (DN_DX'), x, 0,  h);
% Internal forces. Water pressure forces
ElementalMatrix([1, 3],[2, 4]) = -int( DN_DX * N', x, 0,  h);

% Mass conservation
% Volume change
ElementalMatrix([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
% Darcy law
ElementalMatrix([2,4],[2,4]) = dt*int( DN_DX*k*DN_DX', x, 0, h);
%Biot Coefficient
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + AlphaEval*h/12*[1, -1; -1 1];

% Matrix to multiply the old time
% Volume change
ElementalMatrix2([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
%Biot Coefficient
ElementalMatrix2([2,4],[2,4]) = ElementalMatrix2([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix2([2,4],[2,4]) = ElementalMatrix2([2,4],[2,4]) + AlphaEval*h/12*[1, -1; -1 1];

ElementalMatrix = eval(ElementalMatrix);
ElementalMatrix2 = eval(ElementalMatrix2);
% Create the system matrix
SystemMatrix =sparse(2*nNodes, 2*nNodes);
SystemMatrix2 = sparse(2*nNodes, 2*nNodes);

for i = 1:(nNodes-1)
    SystemMatrix( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        SystemMatrix( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + ElementalMatrix;
    SystemMatrix2( 2*(i-1) + [1:4], 2*(i-1) + [1:4]) = ...
        SystemMatrix2( 2*(i-1) + [1:4] , 2*(i-1) + [1:4]) + ElementalMatrix2;
end

% Apply dirichlet conditions
% ZeroWaterPressure
SystemMatrix(2, :) = 0;
SystemMatrix(2,2) = 1;

% Displacement
nn = 2*(nNodes -1)+1;
SystemMatrix(nn, :) = 0;
SystemMatrix(nn,nn) = 1;

% Apply dirichlet to the other matrix
SystemMatrix2(2,:) = 0;
SystemMatrix2(nn,:) = 0;


%Right hand side (water pressure at dirichlet)
FFExt = sparse(2*nNodes, 1);
FFExt(2) = ProblemData.FixedWP;
FFExt(1) = ProblemData.LOAD;


% Solve the linear systems
SolutionOld = zeros(2*nNodes, 1);
SolutionOld(4:2:end) = 1;
for t = 1:nTimeSteps
    Solution = SystemMatrix\(FFExt + SystemMatrix2 * SolutionOld);
    SolutionOld = Solution;
end


NodalWaterPressure = (zeros(nNodes, 1));
NodalDisplacement = (zeros(nNodes, 1));
for i = 1:nNodes
    NodalDisplacement(i) = Solution(2*(i-1)+1);
    NodalWaterPressure(i) = -Solution(2*(i-1)+2);
end



TT =  M * tFinal*k;

ProblemData.NodalDisplacement = NodalDisplacement;
ProblemData.NodalWaterPressure = NodalWaterPressure;
ProblemData.NodalWaterDisplacement = [];
ProblemData.TT = TT;
ProblemData.XX =XX;

ProblemData.OrderU = 1;
ProblemData.OrderPW = 1;
ProblemData.xxPW = [];
ProblemData.xxU = [];
