% Solver for a linear problem

function [X, GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 1, dt, true, AlphaStabM);

[C, K, X, f, fini, nDir] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisIProblem-', ElementType]);

A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii-dt*A;
B = inv(B);


B(nDir,nDir) = eye(length(nDir));
invC = inv(C);
invCfini = (C\fini);
invC(nDir,nDir) = 0;
invCfini(nDir) = 0;

[f] = ComputeForceVector(dt, Nodes, Elements, GPInfo);
X = B*(X + invCfini + (1/nSteps)*invC*f);

for i = 2:nSteps
    [f] = ComputeForceVector(i*dt, Nodes, Elements, GPInfo);
    
    X = B*(X +  (1/nSteps)*invC*f);
end


GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);

PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisIProblem-', ElementType]);

