% Solver for a linear problem

function [X, GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, trash)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, trash, dt, true);

[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisIProblem-', ElementType]);

A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii-dt*A;
B = inv(B);



invCf = (C\f);
invCfini = (C\fini);

X = B*(X + invCfini);
for i = 2:nSteps

    X = B*(X + (1/nSteps)* invCf);
end


GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(GPInfo);

PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisIProblem-', ElementType]);

