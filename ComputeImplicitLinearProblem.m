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
invCf = (C\f);
invCfini = (C\fini);
invCf(nDir) = 0;
invCfini(nDir) = 0;
X = B*(X + invCfini);
for i = 2:nSteps

    X = B*(X + (1/nSteps)* invCf);
end


GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(GPInfo);

PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisIProblem-', ElementType]);

