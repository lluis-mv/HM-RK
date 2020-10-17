% Solver for a linear problem

function [X, GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

[C, K, GPInfo, X, f] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

A = C\(K);

%[l, u] = lu(C); 

ii = eye(3*nNodes, 3*nNodes);
B = ii+dt*A;


%X2 = X;

invCf = (C\f);

for i = 1:nSteps
%     X2 = dt*(u\(l\(K*X2)))+X2;
    X = B*X + (1/nSteps)* invCf;
end
% Linear, 
GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);
GPInfo = FinalizeConstitutiveLaw(GPInfo);


