
function [ProblemData] = ConstructAndSolveProblemUWwP( ProblemData )

nNodes = ProblemData.nNodes;
tFinal = ProblemData.tFinal;
nTimeSteps = ProblemData.nTimeSteps;
StabMethod = ProblemData.StabMethod;
k = ProblemData.k;
M = ProblemData.M;
nu = ProblemData.nu;
QBiot = ProblemData.QBiot;


dt = tFinal/nTimeSteps;
OneOverBiotQ = 1/QBiot;

h = 1/(nNodes-1);

XX = linspace(0,1,nNodes)';


beta = 0.3025;
gamma = 0.6;


alpha1 = 1 / beta / (dt^2);
alpha2 = 1 / beta / dt;
alpha3 = 1/2/beta -1;
alpha4 = gamma / beta/dt;
alpha5 = 1 - gamma/beta;
alpha6 = (1-gamma/beta/2)*dt;

AlphaStab = 0;
if (StabMethod == 1)
    AlphaStab = 0;
elseif ( StabMethod == 2)
    G = M*(1-2*nu)/2/(1-nu);
    AlphaStab = 1 / 2 / G;
elseif (StabMethod == 3)
    AlphaStab = 2/M - 6 * dt * k / h^2;
    AlphaStab = AlphaStab * (0.5 + 0.5*tanh( 2 - 12*dt * k * M/h^2) );
elseif ( StabMethod == 4)
    AlphaStab = 3/M;
elseif ( StabMethod == 5)
    AlphaStab = 2/M - beta/gamma*12 * dt * k / h^2;
    hCrit = sqrt(6*M*beta*dt*k/gamma )
    if (AlphaStab< 0)
        AlphaStab = 0;
    end
%     AlphaStab = 2/M
elseif ( StabMethod == 6)
    AlphaStab = 3/M - 12 * dt * k / h^2;
    if (AlphaStab< 0)
        AlphaStab = 0;
    end
else
    error 'no stabilization method';
end



if ( AlphaStab < 0)
    AlphaStab = 0;
end


rho_mixture = 0.0;
rho_water = 0.0;
porosity = 0.3;

% Define shape functions (to then integrate)
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
SystemMatrix = sparse(3*nNodes, 3*nNodes);
DampingMatrixBig = sparse(3*nNodes, 3*nNodes);
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


FFExt = zeros(3*nNodes, 1);
FFExt(1) = ProblemData.LOAD;
FFExt(3) = ProblemData.FixedWP;
FF = FFExt;

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
    
end

NodalWaterPressure = zeros(nNodes, 1);
NodalDisplacement = NodalWaterPressure;
NodalWaterDisplacement = NodalWaterPressure;

for i = 1:nNodes
    NodalDisplacement(i) = Solution(3*(i-1)+1);
    NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
    NodalWaterPressure(i) = -Solution(3*(i-1)+3);
end

TT =  M * tFinal*k;

ProblemData.NodalDisplacement = NodalDisplacement;
ProblemData.NodalWaterPressure = NodalWaterPressure;
ProblemData.NodalWaterDisplacement = NodalWaterDisplacement;
ProblemData.TT = TT;
ProblemData.XX =XX;


