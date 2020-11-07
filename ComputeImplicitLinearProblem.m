% Solver for a linear problem

function [X, GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, dt, true);

[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);

[C2, K2 ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP,  ElementType, dt, false, 1, 2, 0);

[C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K2);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);

A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii-dt*A;
B = inv(B);

invCf = (C2\f);
invCfini = (C2\fini);

X = B*(X + invCfini);
for i = 2:nSteps

    X = B*(X + (1/nSteps)* invCf);
end


GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(GPInfo);

PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);

