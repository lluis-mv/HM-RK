
% Solver for a linear problem

function [X, GPInfo, normResidual] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, drift)


drift = false;



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, dt, false);

[C, K, X, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);


ii = eye(3*nNodes, 3*nNodes);



[GPInfo] = InitializeConstitutiveLaw(GPInfo);
GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);

if (nargout == 3)
    [f0] = ComputeInternalForces(Elements, GPInfo, X);
end

UnbalancedForces = 0*f0;

for i = 1:nSteps
    
    % Create again C with the appropriate ElastoPlastic stiffness matrix
    [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType,  dt, false);

    [C, K,  ~, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
    
    [C2, ~ ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, dt, false, 1, 2, 0);

    [C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K);
    
    %[C3, ~ ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, dt, false, 1, drift, 0);

    %[C3, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C3, K);
    %C2 = C;
    
    A = C\(K);
    B = ii+dt*A;
    
    invCf = (C2\( (1/nSteps)*f));
    
    
    if ( drift)
        invCf2 = C3\(UnbalancedForces) ;
        invCf = invCf+invCf2;
    end
        
    if ( i == 1)
        invCfini = (C2\fini);
        X = B*X + invCfini+invCf;
    else  
        X = B*X + invCf;
    end
    
    
    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    PostProcessResults(Nodes, Elements, X, GPInfo, dt*i, false, ['ThisProblem-', ElementType]);
    
    if ( drift)
        UnbalancedForces = fini + f*(i/nSteps) + f0 - ComputeInternalForces( Elements, GPInfo, X);
        UnbalancedForces(nDirichlet) = 0;
    end
    
end

if ( drift)
    invCf2 = C3\(UnbalancedForces) ;
    X = X + invCf2;
    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
end


PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);



% Compute the mechanical residual, just to have some fun,....
if (nargout == 3)
    % mechanical part
	finter = ComputeInternalForces( Elements, GPInfo, X);
   
    residual = fini + f*(i/nSteps) + f0 - finter;
    residual(nDirichlet) = 0;
    normResidual = norm(residual);
end
