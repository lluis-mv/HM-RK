
% Solver for a linear problem

function [X, GPInfo, normResidual] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, drift)

if (nargin == 5)
    drift = false;
end


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

[C, K, X, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);


ii = eye(3*nNodes, 3*nNodes);



[GPInfo] = InitializeConstitutiveLaw(GPInfo);
GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);

if (nargout == 3)
    [f0] = ComputeInternalForces(Elements, GPInfo, X);
end

UnbalancedForces = 0*f0;

for i = 1:nSteps
    
    % Create again C with the appropriate ElastoPlastic stiffness matrix
    [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt);

    [C, K,  ~, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
    
    [C2, ~ ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, false, 1, 2, 0);

    [C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K);
    %C2 = C;
    
    A = C\(K);
    B = ii+dt*A;
    
    invCf = (C2\( (1/nSteps)*f));
    
    
    if ( drift)
        invCf2 = C2\(UnbalancedForces) ;
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
    
    if ( drift)
        UnbalancedForces = fini + f*(i/nSteps) + f0 - ComputeInternalForces( Elements, GPInfo, X);
        UnbalancedForces(nDirichlet) = 0;
    end
    
end

if ( drift)
    invCf2 = C2\(UnbalancedForces) ;
    X = X + invCf2;
    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
end



% Compute the mechanical residual, just to have some fun,....
if (nargout == 3)
    % mechanical part
	finter = ComputeInternalForces( Elements, GPInfo, X);
   
    residual = fini + f*(i/nSteps) + f0 - finter;
    residual(nDirichlet) = 0;
    normResidual = norm(residual);
end
