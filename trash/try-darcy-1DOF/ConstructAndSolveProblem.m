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



% Define shape functions (to then integrate)
syms x positive
chi = 2*x/h-1;
if ( ProblemData.OrderU == 1)
    N = [(h-x)/h; x/h];
elseif (ProblemData.OrderU == 2)
    N = [1/2*chi*(chi-1);  (1-chi^2); 1/2*chi*(chi+1)];
elseif ( ProblemData.OrderU == 3)
    N = [ -9/16*(chi+1/3)*(chi-1/3)*(chi-1); %N1
    27/16*(chi+1)*(chi-1/3)*(chi-1); %N2
    -27/16*(chi+1)*(chi+1/3)*(chi-1); %N3
    9/16*(chi+1)*(chi+1/3)*(chi-1/3) %N4
    ];
elseif (ProblemData.OrderU == 4)
               N = [ (2*chi*(chi - 1)*(chi - 1/2)*(chi + 1/2))/3;
    -(8*chi*(chi - 1)*(chi + 1)*(chi - 1/2))/3;
 4*(chi - 1)*(chi + 1)*(chi - 1/2)*(chi + 1/2);
    -(8*chi*(chi - 1)*(chi + 1)*(chi + 1/2))/3;
   (2*chi*(chi + 1)*(chi - 1/2)*(chi + 1/2))/3;
   ];
end

DN_DX = diff(N, x);

nDofs = ProblemData.OrderU+1;
% Define the elemental system
ElementalMatrix = sym(zeros(nDofs,nDofs));



% Internal forces. Effective stress forces
ElementalMatrix = -int(DN_DX * (DN_DX'), x, 0,  h);


ElementalMass = int(N * N', x, 0,  h);
% Lump the matrix
EE = 0 * ElementalMass;
for i = 1:nDofs
    EE(i,i) = sum(ElementalMass(:,i));
end
% ElementalMass = EE;

ElementalMatrix = eval(ElementalMatrix);
% Create the system matrix
nSystemSize = nNodes + (nDofs-2)*(nNodes-1);

SystemMatrix =sparse(nSystemSize,nSystemSize);
MassMatrix = SystemMatrix;
FFExt = sparse(nSystemSize, 1);

clear x;
clear N;

syms x real
for i = 1:(nNodes-1)
    ind = (i-1)*(nDofs-1) + [1:nDofs];
    SystemMatrix(ind,ind) = SystemMatrix( ind,ind) + ElementalMatrix;
    MassMatrix(ind,ind) = MassMatrix( ind,ind) + ElementalMass;

    % TryToEvaluateSourceTerm
    if (ProblemData.tFinal> 1)
        xIni = XX(i); xFin = XX(i+1);
        chi = 2*(x-xIni)/h-1;
        if ( ProblemData.OrderU == 1)
        chi = x - xIni;
        N = [(h-chi)/h; chi/h];
    elseif (ProblemData.OrderU == 2)
        N = [1/2*chi*(chi-1);  (1-chi^2); 1/2*chi*(chi+1)];
    elseif ( ProblemData.OrderU == 3)
        N = [ -9/16*(chi+1/3)*(chi-1/3)*(chi-1); %N1
        27/16*(chi+1)*(chi-1/3)*(chi-1); %N2
        -27/16*(chi+1)*(chi+1/3)*(chi-1); %N3
        9/16*(chi+1)*(chi+1/3)*(chi-1/3) %N4
        ];
    elseif (ProblemData.OrderU == 4)
               N = [ (2*chi*(chi - 1)*(chi - 1/2)*(chi + 1/2))/3;
    -(8*chi*(chi - 1)*(chi + 1)*(chi - 1/2))/3;
 4*(chi - 1)*(chi + 1)*(chi - 1/2)*(chi + 1/2);
    -(8*chi*(chi - 1)*(chi + 1)*(chi + 1/2))/3;
   (2*chi*(chi + 1)*(chi - 1/2)*(chi + 1/2))/3;
   ];
    end
    FF = int( N*source(x), xIni, xFin);
    FF = eval(FF);
    FFExt(ind) = FFExt(ind) + FF;
    end
end



MMatrix = MassMatrix;
KMatrix = SystemMatrix;
nn = 1;
FFExt(1) = 0;
if ( tFinal > 1)
    FFExt(end)=1;
    
    nn = [1,nSystemSize];
end
% Displacement


SystemMatrix(nn, :) = 0;
SystemMatrix(nn,nn) = eye(length(nn));

MassMatrix(nn,:)=0;

%Right hand side (water pressure at dirichlet)


if ( dt < 1)
%     SolutionOld = ones(nSystemSize,1);
%     SolutionOld(1) = 0;
%     alfa = 1;
%     SystemMatrix(nn,nn) = -1./dt/alfa * SystemMatrix(nn,nn);
% %     FFExt(1) =1;
%     
%     for i = 1:ProblemData.nTimeSteps
%         Solution = (MassMatrix - dt *alfa * SystemMatrix) \ (FFExt + (MassMatrix +(1-alfa)*dt*SystemMatrix)  * SolutionOld);
%         SolutionOld = Solution;
%     end
% 
%     SolutionOld = ones(nSystemSize,1);
%     SolutionOld(1) = 0;
    
    Sys = MMatrix \ KMatrix;
    Sys(1,:) = 0;
    Sys(1,1) = 1;
%     opts = odeset('RelTol',1e-13,'AbsTol',1e-14);
%     [t, y] = ode45(@(t,y) Sys*y, linspace(0, ProblemData.tFinal), SolutionOld, opts);
%     Solution = y(end,:)';

    % Now Try To Integrate 
    SolutionOld = ones(nSystemSize,1);
    SolutionOld(1) = 0;
    [vectors, values] = eig(full(Sys));
    Constants = vectors\SolutionOld;
    
    tt= linspace(-5,1,20)
    tt = 10.^tt;
    tt = ProblemData.tFinal
    YY = [];
    for t = tt
        ySol = zeros(nSystemSize,1);
        for n = 1:nSystemSize
            ySol = ySol + Constants(n)*exp(values(n,n)*t) * vectors(:,n);
        end
        YY = [YY, ySol];
    end
    figure(200)
%     semilogx(tt, YY);
    pause(0.01)
    Solution = ySol;
else
    Solution = SystemMatrix\(FFExt);
end





ProblemData.NodalDisplacement = full(Solution)';
ProblemData.xxU = linspace(0,1,length(Solution));
ProblemData.XX =XX;


