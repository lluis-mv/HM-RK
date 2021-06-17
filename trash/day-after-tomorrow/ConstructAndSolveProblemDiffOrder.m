function ProblemData = ConstructAndSolveProblemDiffOrder( ProblemData)


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
Np = [(h-x)/h; x/h];
DN_DXP = diff(Np, x);


chi = 2*x/h-1;
Nu = [1/2*chi*(chi-1); 1/2*chi*(chi+1); (1-chi^2)];
DN_DXU = diff(Nu, x);

% Define the elemental system
ElementalMatrix = sym(zeros(5,5));
ElementalMatrix2 = sym(zeros(5,5));

indexu = [1,2,3];
indexp = [4,5];

% Internal forces. Effective stress forces
ElementalMatrix(indexu, indexu) = int( M *  DN_DXU * (DN_DXU'), x, 0,  h);
% Internal forces. Water pressure forces
ElementalMatrix(indexu, indexp) = -int( DN_DXU * Np', x, 0,  h);

% Mass conservation
% Volume change
ElementalMatrix(indexp,indexu) = int( Np*DN_DXU', x, 0, h);
% Darcy law
ElementalMatrix(indexp,indexp) = dt*int( DN_DXP*k*DN_DXP', x, 0, h);
%Biot Coefficient
ElementalMatrix(indexp,indexp) = ElementalMatrix(indexp,indexp) + (1/QBiot)*int( Np*Np', x, 0, h);



% Matrix to multiply the old time
% Volume change
ElementalMatrix2(indexp,indexu) = int( Np*DN_DXU', x, 0, h);
%Biot Coefficient
ElementalMatrix2(indexp,indexp) = ElementalMatrix2(indexp,indexp) + (1/QBiot)*int( Np*Np', x, 0, h);


ElementalMatrix = eval(ElementalMatrix);
ElementalMatrix2 = eval(ElementalMatrix2);

% Reshape The Matrices to have a similar form than before
index = [1,4,2,5,3];
ElementalMatrix = ElementalMatrix(index,index);
ElementalMatrix2 = ElementalMatrix2(index,index);


% Create the system matrix
nElements = nNodes-1;
nSystemSize = 2*nNodes + nElements;

SystemMatrix =sparse(nSystemSize, nSystemSize);
SystemMatrix2 = sparse(nSystemSize,nSystemSize);

last = 2*nNodes+1;
for i = 1:nElements
    index = [2*(i-1)+[1:4], last];
    SystemMatrix(index, index) = ...
        SystemMatrix( index,index) + ElementalMatrix;
    SystemMatrix2( index,index) = ...
        SystemMatrix2( index,index) + ElementalMatrix2;
    last = last+1;
end

% Apply dirichlet conditions
% ZeroWaterPressure
nFixWaterPressure = 2;
SystemMatrix(nFixWaterPressure,:) = 0;
SystemMatrix(nFixWaterPressure,nFixWaterPressure) = 1;

% Displacement
nn = 2*(nNodes -1)+1;
SystemMatrix(nn, :) = 0;
SystemMatrix(nn,nn) = 1;

% Apply dirichlet to the other matrix
SystemMatrix2(nFixWaterPressure,:) = 0;
SystemMatrix2(nn,:) = 0;


%Right hand side (water pressure at dirichlet)
FFExt = sparse(nSystemSize, 1);
FFExt(nFixWaterPressure) = ProblemData.FixedWP;
FFExt(1) = ProblemData.LOAD;


% Solve the linear systems
SolutionOld = zeros(nSystemSize, 1);
SolutionOld(4:2:2*nNodes) = 1;
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


last = 2*nNodes+1;
j = 1;
for i = 1:nElements
    xxU(j) = XX(i);
    solU(j) = NodalDisplacement(i);
    j = j+1;
    xxU(j) = 0.5*(XX(i)+XX(i+1));
    solU(j) = Solution(last);
    last = last+1;
    j = j+1;
end
xxU = [xxU, XX(end)];
solU = [solU, NodalDisplacement(end)];

ProblemData.NodalDisplacement = solU;
ProblemData.xxU = xxU;

ProblemData.OrderU = 2;
ProblemData.OrderPW = 1;
ProblemData.xxPW = [];
