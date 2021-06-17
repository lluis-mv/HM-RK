function [ProblemData] = SolveAndComputeNorms( ProblemData)


if ( ProblemData.PrimalForm)
    if ( ProblemData.StabMethod < 6)
        ProblemData = ConstructAndSolveProblem( ProblemData );
    elseif ( ProblemData.StabMethod ==6)
        ProblemData = ConstructAndSolveProblemDiffOrder( ProblemData);
    elseif ( ProblemData.StabMethod == 7)
        ProblemData = ConstructAndSolveProblemDiffOrderP3P1( ProblemData);
    elseif ( ProblemData.StabMethod == 8)
        ProblemData = ConstructAndSolveProblemDiffOrderP3P2( ProblemData);
    else
        error
    end
else
    ProblemData = ConstructAndSolveProblemUWwP( ProblemData );
end

nGP = 3; % i have to ask someone that knows things how many gauss points are necessary to compute the error
XX = ProblemData.XX;
TT = ProblemData.TT;
NodalWaterPressure = ProblemData.NodalWaterPressure;
NodalDisplacement = ProblemData.NodalDisplacement;

Draw = false;

[xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
L2 = EvaluateL2( XX, pwGP, NodalWaterPressure, nGP);
if ( Draw)
    figure(103)
    xPlot = [];
    yPlot = [];
    for a = 1:size(xGP, 1)
        
        xPlot = [xPlot, xGP(a,:)];
        yPlot = [yPlot, pwGP(a,:)];
    end
    plot(xPlot, yPlot, 'r')
    hold on
end

pw = EvalConsolidation(XX, TT);
LInf = max( abs(NodalWaterPressure - pw));
if ( Draw)
    plot(XX, pw, 'k')
    plot(XX, NodalWaterPressure, 'b')
    hold off
end

[xGP, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP);
L2DISPL = EvaluateL2( XX, uuGP, NodalDisplacement, nGP);

if ( Draw)
    figure(104)
    xPlot = [];
    yPlot = [];
    for a = 1:size(xGP, 1)
        xPlot = [xPlot, xGP(a,:)];
        yPlot = [yPlot, uuGP(a,:)];
    end
    plot(xPlot, yPlot, 'r')
    hold on
end

uu = EvalConsolidationDISPL(XX, TT);
LInfDISPL = max( abs(NodalDisplacement - uu));
if ( Draw)
    plot(XX, uu, 'k')
    plot(XX, NodalDisplacement, 'b')
    hold off
    if ( isempty(ProblemData.NodalWaterDisplacement) ==false )
        figure(105)
        plot(XX, ProblemData.NodalWaterDisplacement, 'b')
    end
    pause(1)
end

ProblemData.LInfDISP = LInfDISPL;
ProblemData.L2DISP = L2DISPL;

ProblemData.LInf = LInf;
ProblemData.L2 = L2;




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
ProblemData.TT = TT;
ProblemData.XX =XX;



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
    if (AlphaStab< 0)
        AlphaStab = 0;
    end
else
    AlphaStab = 3/M - 12 * dt * k / h^2;
    if (AlphaStab< 0)
        AlphaStab = 0;
    end
end



if ( AlphaStab < -1e-8)
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



function pw = EvalConsolidation(XX, TT)
pw = 0*XX;
nTermsZero = 0;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:5000
            aux = pi/2*(2*m+1);
            term = 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT);
            if ( abs(term) < 1e-12)
                nTermsZero = nTermsZero + 1;
            else 
                nTermsZero = 0;
            end
            pw(i,j) = pw(i,j)+term ;
            if ( nTermsZero > 20)
                break;
            end
        end
    end
end

function uu = EvalConsolidationDISPL(XX, TT)
uu = -(1-XX);
nTermsZero = 0;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:5000
            z=XX(i,j);
            term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
                   if ( abs(term) < 1e-12)
                nTermsZero = nTermsZero + 1;
            else 
                nTermsZero = 0;
            end
            uu(i,j) = uu(i,j)+term;
            if ( nTermsZero > 20)
                break;
            end
        end
    end
end

function [xPG, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP)

if ( nGP == 1)
    xPG = [];
    for i = 1:length(XX)-1
        xPG = [xPG; mean(XX(i:i+1))];
    end
elseif (nGP == 2)
    xPG=[];
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
elseif (nGP == 3)
    xPG=[];
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
end

uuGP = EvalConsolidationDISPL(xPG, TT);

function [xPG, pwPG] = EvalConsolidationGP(XX, TT, nGP)

if ( nGP == 1)
    xPG = [];
    for i = 1:length(XX)-1
        xPG = [xPG; mean(XX(i:i+1))];
    end
elseif (nGP == 2)
    xPG=[];
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
elseif (nGP == 3)
    xPG=[];
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    for i = 1:length(XX)-1
        xx = N*[XX(i); XX(i+1)];
        xPG = [xPG; xx'];
    end
end


pwPG = EvalConsolidation(xPG, TT);

function subError = EvaluateL2( XX, pwGP, WPNodes, nGP)

subError = 0;

if ( nGP == 1)
    L = XX(2)-XX(1);
    w = L/2*2;
    for i = 1:length(XX)-1
        N = [0.5, 0.5];
        nodal = N*(WPNodes(i:i+1));
        subError = subError + (nodal-pwGP(i))^2*w;
    end
elseif (nGP == 2)
    N = [0.5+1/2/sqrt(3), 0.5-1/2/sqrt(3);
        0.5-1/2/sqrt(3), 0.5+1/2/sqrt(3)];
    L = XX(2)-XX(1);
    w = L/2;
    for i = 1:length(XX)-1
        nodal = N * [WPNodes(i:i+1)];
        pwLocal = pwGP(i,:);
        for a = 1:2
            subError = subError + (nodal(a) -pwLocal(a))^2 * w;
        end
    end
elseif (nGP == 3)
    N = [0.5+1/2*sqrt(3/5), 0.5, 0.5-1/2*sqrt(3/5);
        0.5-1/2*sqrt(3/5), 0.5, 0.5+1/2*sqrt(3/5)]';
    L = XX(2)-XX(1);
    w = L/2*[5/9, 8/9,5/9];
    for i = 1:length(XX)-1
        nodal = N * [WPNodes(i:i+1)];
        pwLocal = pwGP(i,:);
        for a = 1:3
            subError = subError + (nodal(a) -pwLocal(a))^2 * w(a);
        end
    end
end

subError = sqrt(subError);
