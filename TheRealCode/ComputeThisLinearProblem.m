% Solver for a linear problem

function [X, GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, false);

[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);

[C2, K2 ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, false, 1, 2, 0);

[C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K2);


A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii+dt*A;

invCf = (C2\f);
invCfini = (C2\fini);

X = B*X+invCfini;
hola = 1;
for i = 2:nSteps
    X = B*X + (1/nSteps)* invCf;
end



% Linear, 
GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);
GPInfo = FinalizeConstitutiveLaw(GPInfo);


