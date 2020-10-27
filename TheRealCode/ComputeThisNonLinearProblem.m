% Solver for a linear problem

function [X, GPInfo] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);


ii = eye(3*nNodes, 3*nNodes);

GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);

GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);
[GPInfo] = InitializeConstitutiveLaw(GPInfo);

for i = 1:nSteps
    
    % Create again C with the appropriate ElastoPlastic stiffness matrix
    [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

    [C, K,  ~, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
    
    [C2, ~ ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, false, 1, 2, 0);

    [C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K);

    
    A = C\(K);
    B = ii+dt*A;
    
    invCf = (C\f);
    
    if ( i == 1)
        invCfini = (C2\fini);
        X = B*X + invCfini;
    else  
        X = B*X + (1/nSteps)* invCf;
    end

    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    
end
