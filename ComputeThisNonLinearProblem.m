
% Solver for a linear problem

function [X, GPInfo, normResidual] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false);

[C, K, X, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);



if ( RKMethod)
    [a,b] = GetRungeKutta(RKMethod);
end


[GPInfo] = InitializeConstitutiveLaw(GPInfo);
GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);

if (nargout == 3)
    [f0] = ComputeInternalForces(Elements, GPInfo, X);
end


for t = 1:nSteps
    
    for i = 1:length(b)
        XStep = X;
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
        GPInfo = EvaluateConstitutiveLaw(GPInfo, XStep, Elements, false);
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod,  dt, false);

        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C, K);
    
        A = C\K;
        invCf = (C\f);
        k(:,i) = A*XStep + (1/(dt*nSteps))*invCf;
        if ( t == 1)
            k(:,i) = k(:,i) + (1/dt)* (C\fini);
        end
    end
    XNew  = X; 
   for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end     
  
    X = XNew;
    
    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    %PostProcessResults(Nodes, Elements, X, GPInfo, dt*i, false, ['ThisProblem-', ElementType]);
    
%     if ( drift)
%         UnbalancedForces = fini + f*(i/nSteps) + f0 - ComputeInternalForces( Elements, GPInfo, X);
%         UnbalancedForces(nDirichlet) = 0;
%     end
    
end



PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);



% Compute the mechanical residual, just to have some fun,....
if (nargout == 3)
    % mechanical part
	finter = ComputeInternalForces( Elements, GPInfo, X);
   
    residual = fini + f + f0 - finter;
    residual(nDirichlet) = 0;
    normResidual = norm(residual);
end
