function []=main5()

clc; clear all; close all;

if (false)
    nStab = 5;
    i = 1;
    
    % Effect of Poisson Ratio
    nNodes = 20;
    TEND = 0.0001;
    
    % Figure 1. Effect of nu
    NU = [0.00:0.1:0.40 0.45 0.48, 0.49, 0.499, 0.4999];
    NU = [NU, 1/3];
    NU = sort(NU)
    Draw = false;
    for nu = NU
        for StabMethod = [1:nStab]
            if ( Draw)
                figure(3)
                subplot(3,2,StabMethod)
                figure(4)
                subplot(3,2,StabMethod)
            end
            [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblem(20, 0.0001, 1, StabMethod, Draw, nu);
            L2Matrix(i,StabMethod) = L2;
            LInfMatrix(i,StabMethod) = LInf;
            L2MatrixDISPL(i,StabMethod) = L2DISPL;
            LInfMatrixDISPL(i,StabMethod) = LInfDISPL;
            if ( Draw)
                figure(3)
                plot(linspace(0,1,nNodes), WP)
                hold off
                figure(4)
                plot(linspace(0,1,nNodes), U)
                hold off
            end
        end
        pause(0.0000001)
        i = i+1
    end
    
    XBAR = 0.5-NU;
    figure(1)
    for i = 1:nStab
        loglog(XBAR, L2Matrix(:,i), '*-.')
        hold on
    end
    hold on
    ylabel('$\| p-p_h\|_{L2}$', 'interpreter', 'latex')
    xlabel('$0.5-\nu$', 'interpreter', 'latex')
    
    
    figure(2)
    for i = 1:nStab
        loglog(XBAR, LInfMatrix(:,i), '*-.')
        hold on
    end
    hold on
    ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
    xlabel('$0.5-\nu$', 'interpreter', 'latex')
    
    
    
    figure(11)
    for i = 1:nStab
        loglog(XBAR, L2MatrixDISPL(:,i), '*-.')
        hold on
    end
    hold on
    ylabel('$\| u-u_h\|_{L2}$', 'interpreter', 'latex')
    xlabel('$0.5-\nu$', 'interpreter', 'latex')
    
    
    figure(12)
    for i = 1:nStab
        loglog(XBAR, LInfMatrixDISPL(:,i), '*-.')
        hold on
    end
    hold on
    ylabel('$\|u(x_i) - u^h_i\|$ (kPa)','interpreter', 'latex')
    xlabel('$0.5-\nu$', 'interpreter', 'latex')
    
    
    for i = [1,2,11,12]
        figure(i)
        ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', 'LMV-2','location', 'best');
        set(ll, 'interpreter', 'latex');
        print(['IMAGE-POISSON', num2str(i)], '-dpng')
    end
    
    
    
end


if ( false)
    nStab = 5;
    
    i=1; L2Matrix =[]; L2MatrixDISPL=[]; LInfMatrix=[]; LInfMatrixDISPL=[];
    nNodes = 20;
    TEND = 0.0001;
    
    % Figure 1. Effect of nu
    TEND = [-8:0.5:-1];
    TEND = 10.^TEND
    Draw = false;
    for te = TEND
        for StabMethod = [1:nStab]
            if ( Draw)
                figure(3)
                subplot(3,2,StabMethod)
                figure(4)
                subplot(3,2,StabMethod)
            end
            [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblem(nNodes, te, 1, StabMethod);
            L2Matrix(i,StabMethod) = L2;
            LInfMatrix(i,StabMethod) = LInf;
            L2MatrixDISPL(i,StabMethod) = L2DISPL;
            LInfMatrixDISPL(i,StabMethod) = LInfDISPL;
            if ( Draw)
                figure(3)
                plot(linspace(0,1,nNodes), WP)
                hold off
                figure(4)
                plot(linspace(0,1,nNodes), U)
                hold off
            end
        end
        pause(0.0000001)
        i = i+1
    end
    
    XBAR = TEND;
    figure(1)
    for i = 1:nStab
        loglog(XBAR, L2Matrix(:,i), '*-.')
        hold on
    end
    hold on
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 dt k/ h^2 $)', 'LMV-2','location', 'best');
    set(ll, 'interpreter', 'latex');
    ylabel('$\| p-p_h\|_{L2} (kPa)$', 'interpreter', 'latex')
    xlabel('$\Delta t$', 'interpreter', 'latex')
    
    
    figure(2)
    for i = 1:nStab
        loglog(XBAR, LInfMatrix(:,i), '*-.')
        hold on
    end
    hold on
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', 'LMV-2','location', 'best');
    set(ll, 'interpreter', 'latex');
    ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
    xlabel('$\Delta t$', 'interpreter', 'latex')
    
    
    
    figure(11)
    for i = 1:nStab
        loglog(XBAR, L2MatrixDISPL(:,i), '*-.')
        hold on
    end
    hold on
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 dt k/ h^2 $)', 'LMV-2','location', 'best');
    set(ll, 'interpreter', 'latex');
    ylabel('$\| u-u_h\|_{L2}$ (m)', 'interpreter', 'latex')
    xlabel('$\Delta t$', 'interpreter', 'latex')
    
    
    figure(12)
    for i = 1:nStab
        loglog(XBAR, LInfMatrixDISPL(:,i), '*-.')
        hold on
    end
    hold on
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', 'LMV-2','location', 'best');
    set(ll, 'interpreter', 'latex');
    ylabel('$\|u(x_i) - u^h_i\|$ (m)','interpreter', 'latex')
    xlabel('$\Delta t$', 'interpreter', 'latex')
    
    
end



nStab = 7;
i = 1;

% Effect of Poisson Ratio
nNodes = 20;
TEND = 0.0001;

% Figure 1. Effect of nu
nNodes = [4,6,8,10,15,20,38,39,40,41,42,60,80,100]%, 150, 200]
nNodes = [4:2:80];
nNodes = [4, 6, 8, 10, 20, 40, 80, 200];
size(nNodes)
Draw = false;

L2Matrix = zeros(length(nNodes), nStab);
LInfMatrix = zeros(length(nNodes), nStab);
L2MatrixDISPL = zeros(length(nNodes), nStab);
LInfMatrixDISPL = zeros(length(nNodes), nStab);

for nodes = nNodes
    for StabMethod = [1, 5, 6, 7]
        if ( Draw)
            figure(3)
            subplot(3,2,StabMethod)
            figure(4)
            subplot(3,2,StabMethod)
        end
        if (StabMethod < 6)
            [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblem(nodes, TEND, 1, StabMethod);
        else
            [WP, U, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblemUWwP(nodes, TEND, 1, StabMethod);
        end
        L2Matrix(i,StabMethod) = L2;
        LInfMatrix(i,StabMethod) = LInf;
        L2MatrixDISPL(i,StabMethod) = L2DISPL;
        LInfMatrixDISPL(i,StabMethod) = LInfDISPL;
        if ( Draw)
            figure(3)
            plot(linspace(0,1,nodes), WP)
            hold off
            figure(4)
            plot(linspace(0,1,nodes), U)
            hold off
        end
    end
    pause(0.0000001)
    i = i+1
end

XBAR = 1./(nNodes+1);
figure(1)
for i = 1:nStab
    loglog(XBAR, L2Matrix(:,i), '*-.')
    hold on
end
hold on
ylabel('$\| p-p_h\|_{L2}$', 'interpreter', 'latex')
xlabel('$h_e$', 'interpreter', 'latex')


figure(2)
for i = 1:nStab
    loglog(XBAR, LInfMatrix(:,i), '*-.')
    hold on
end
hold on
ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
xlabel('$h_e$', 'interpreter', 'latex')



figure(11)
for i = 1:nStab
    loglog(XBAR, L2MatrixDISPL(:,i), '*-.')
    hold on
end
hold on
ylabel('$\| u-u_h\|_{L2}$', 'interpreter', 'latex')
xlabel('$h_e$', 'interpreter', 'latex')


figure(12)
for i = 1:nStab
    loglog(XBAR, LInfMatrixDISPL(:,i), '*-.')
    hold on
end
hold on
ylabel('$\|u(x_i) - u^h_i\|$ (kPa)','interpreter', 'latex')
xlabel('$h_e$', 'interpreter', 'latex')


for i = [1,2,11,12]
    figure(i)
    plot(sqrt(TEND*6)*[1,1],   [1e-2 1e0], 'k:')
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', ...
        '$\mathbf{u}-\mathbf{w}-p_w$','$\mathbf{u}-\mathbf{w}-p_w$ (NonStab)', 'location', 'best');
    set(ll, 'interpreter', 'latex');
    print(['IMAGE-DX', num2str(i)], '-dpng')
end















function [NodalWaterPressure, NodalDisplacement, L2, LInf, L2DISPL, LInfDISPL] = ConstructAndSolveProblem( nNodes, tFinal, nTimeSteps, StabMethod, Draw,  nu)

nGP = 3;

syms AlphaStab real

dt = tFinal/nTimeSteps;

QBiot = 1e67;

M =  1;
if ( nargin < 5)
    Draw = false;
end
if ( nargin < 6)
    nu = 1/3;
end

k = 1;

h = 1/(nNodes-1);
G = M*(1-2*nu)/2/(1-nu);
XX = linspace(0,1,nNodes)';

% Define shape functions (to then integrate)
syms x positive
N = [(h-x)/h; x/h];
DN_DX = diff(N, x);

% Define the elemental system
ElementalMatrix = sym(zeros(4,4));
ElementalMatrix2 = sym(zeros(4,4));

% Internal forces. Effective stress forces
ElementalMatrix([1, 3],[1, 3]) = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% Internal forces. Water pressure forces
ElementalMatrix([1, 3],[2, 4]) = -int( DN_DX * N', x, 0,  h);

% Mass conservation
% Volume change
ElementalMatrix([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
% Darcy law
ElementalMatrix([2,4],[2,4]) = dt*int( DN_DX*k*DN_DX', x, 0, h);
%Biot Coefficient
%ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix([2,4],[2,4]) = ElementalMatrix([2,4],[2,4]) + AlphaStab*h/12*[1, -1; -1 1];

% Matrix to multiply the old time
% Volume change
ElementalMatrix2([2,4],[1,3]) = int( N*DN_DX', x, 0, h);
%Biot Coefficient
%ElementalMatrix2([2,4],[2,4]) = ElementalMatrix2([2,4],[2,4]) + (1/QBiot)*int( N*N', x, 0, h);
%Stabilization
ElementalMatrix2([2,4],[2,4]) = ElementalMatrix2([2,4],[2,4]) + AlphaStab*h/12*[1, -1; -1 1];


% Create the system matrix
SystemMatrix = sym( zeros(2*nNodes, 2*nNodes) );
SystemMatrix2 = sym( zeros(2*nNodes, 2*nNodes) );

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
FFExt = sym(zeros(2*nNodes, 1));
FFExt(2) = 0;
FFExt(1) = -1;

if (StabMethod == 1)
    AlphaEval = 0;
    
elseif ( StabMethod == 2)
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
SystemMatrix = subs(SystemMatrix, AlphaStab, AlphaEval);
SystemMatrix = eval(SystemMatrix);

SystemMatrix2 = subs(SystemMatrix2, AlphaStab, AlphaEval);
SystemMatrix2 = eval(SystemMatrix2);

FFExt = eval(FFExt);


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


% Compute Error Norms
TT =  M * tFinal*k;


[xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
L2 = EvaluateL2( XX, pwGP, NodalWaterPressure, nGP);
if ( Draw)
    xPlot = [];
    yPlot = [];
    for a = 1:size(xGP, 1)
        figure(3)
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
end



[xGP, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP);
L2DISPL = EvaluateL2( XX, uuGP, NodalDisplacement, nGP);

if ( Draw)
    figure(4)
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
    hold on
end

function [NodalWaterPressure, NodalDisplacement, L2, LInf, L2DISPL, LInfDISPL]  = ConstructAndSolveProblemUWwP(nNodes, tFinal, nTimeSteps, StabMethod, Draw, nu)

nGP = 3;




dt = tFinal/nTimeSteps;
QBiot = 1e7;
OneOverBiotQ = 1/QBiot;

beta = 0.3025;
gamma = 0.6;

alpha1 = 1 / beta / (dt^2);
alpha2 = 1 / beta / dt;
alpha3 = 1/2/beta -1;
alpha4 = gamma / beta/dt;
alpha5 = 1 - gamma/beta;
alpha6 = (1-gamma/beta/2)*dt;


M =  1;
if ( nargin < 5)
    Draw = false;
end
if ( nargin < 6)
    nu = 1/3;
end

k = 1;

h = 1/(nNodes-1);
G = M*(1-2*nu)/2/(1-nu);
XX = linspace(0,1,nNodes)';


if ( StabMethod == 6)
    AlphaStab = 2/M;
    %AlphaStab = AlphaStab * (0.5 + 0.5*tanh( 2 - 12*dt * k * M/h^2) );
    if (AlphaStab < 0)
        AlphaStab = 0;
    end
else 
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
SystemMatrix = ( zeros(3*nNodes, 3*nNodes) );
DampingMatrixBig = ( zeros(3*nNodes, 3*nNodes) );
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

FFExt = zeros(3*nNodes, 1);

FFExt(3) = 1;
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
    
    if ( Draw)
        for i = 1:nNodes
            NodalDisplacement(i) = Solution(3*(i-1)+1);
            NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
            NodalWaterPressure(i) = 1-Solution(3*(i-1)+3);
        end
        TT =  M *k*step/nTimeSteps * tFinal;
        figure(1);
        subplot(3,1,1)
        uu = EvalConsolidationDISPL(XX, TT);
        plot(XX, uu, 'k')
        hold on
        plot(XX, NodalDisplacement)
        hold off
        subplot(3,1,2)
        plot(XX, NodalWaterDisplacement)
        subplot(3,1,3)
        
        pw = EvalConsolidation(XX, TT);
        plot(XX, pw, 'k')
        hold on
        plot(XX, NodalWaterPressure)
        hold off
    end
end

for i = 1:nNodes
    NodalDisplacement(i) = Solution(3*(i-1)+1);
    NodalWaterDisplacement(i) = Solution(3*(i-1)+2);
    NodalWaterPressure(i) = 1-Solution(3*(i-1)+3);
end

% Compute Error Norms
TT =  M * tFinal*k;


[xGP, pwGP] = EvalConsolidationGP(XX, TT, nGP);
L2 = EvaluateL2( XX, pwGP, NodalWaterPressure, nGP);
if ( Draw)
    xPlot = [];
    yPlot = [];
    for a = 1:size(xGP, 1)
        figure(3)
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
end



[xGP, uuGP] = EvalConsolidationGPDISPL(XX, TT, nGP);
L2DISPL = EvaluateL2( XX, uuGP, NodalDisplacement, nGP);

if ( Draw)
    figure(4)
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
    hold on
end






function pw = EvalConsolidation(XX, TT)
pw = 0*XX;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:100
            aux = pi/2*(2*m+1);
            term = 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT);
            pw(i,j) = pw(i,j)+term ;
        end
    end
end

function uu = EvalConsolidationDISPL(XX, TT)
uu = -(1-XX);
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:100
            z=XX(i,j);
            term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
            uu(i,j) = uu(i,j)+term;
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
