% Solver for a linear problem

function [X, GPInfo] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

[C, K, GPInfo,  X, f] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);


ii = eye(3*nNodes, 3*nNodes);

GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);

for i = 1:nSteps
    
    % Create again C with the appropriate ElastoPlastic stiffness matrix
    [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

    [C, K, GPInfo,  ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
    
    A = C\(K);
    B = ii+dt*A;
    
    invCf = (C\f);
    
    X = B*X + (1/nSteps)* invCf;

    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    
end
