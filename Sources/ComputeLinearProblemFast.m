% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeLinearProblemFast(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

if (nargout == 3)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

if (nargin ==7)
    AlphaStabM = 1;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

fini(nDirichlet) = 0;
A = C\(K);
A(nDirichlet,:)=0;

[~, ~, AllZero] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP);
if ( AllZero)
    invC = 0*C;
else
    invC = inv(C);
end

invCfini = (C\fini);
invCfini(nDirichlet)=0;


if ( RKMethod)
    [a,b, c] = GetRungeKutta(RKMethod);
end
k = zeros(  3*nNodes, length(b));




PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);
if ( DoSomePostProcess )
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end

AmplificationMatrix = ComputeAmplificationMatrix(a, b, dt, A);

for loadStep = 1:nSteps

    X = AmplificationMatrix*X;
    
    if (any(isnan(X)))
        X = nan*X;
        return;
    end
    if (any(isinf(X)))
        X = nan*X;
        return;
    end
    if ( DoSomePostProcess )
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
end




GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false, dt);
GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);


function M = ComputeAmplificationMatrix(a, b, dt, A)
M = (eye(size(A)));
for i = 1:length(a)
    kk = M;
    for j = 1:i-1
        kk = kk + dt*a(i,j)*ki(:,:,j);
    end
    ki(:,:,i) = A * kk;
end

for i = 1:length(a)
    M = M + dt * b(i)*ki(:,:,i);
end

M = sparse(M);
