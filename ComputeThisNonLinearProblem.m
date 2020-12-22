
% Solver for a linear problem

function [X, GPInfo, normResidual, ThisInfo] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM, drift)

if (nargin == 8)
    drift = false;
end
if (nargout == 4)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

RKMethodLaw = RKMethod;
if ( isfield( CP, 'RKMethodLaw'))
    RKMethodLaw = CP.RKMethodLaw;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPInfo] = InitializeConstitutiveLaw(GPInfo);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[~,~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);


if ( any([GPInfo.VonMises] == true) )
    addpath('../ModifiedCamClay/vonMises/')
elseif ( any([GPInfo.MCC] == true) )
    addpath('../ModifiedCamClay/')
end

GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false, RKMethodLaw);
f0 = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
f0(nDirichlet) = 0;

if ( RKMethod)
    [a,b,c] = GetRungeKutta(RKMethod);
end
k = zeros(  3*nNodes, length(b));




PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);
if ( DoSomePostProcess ) 
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end

for loadStep = 1:nSteps
    
    for i = 1:length(b)
        XStep = X;
        t = (loadStep-1)*dt + c(i)*dt;
        [f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
        f(nDirichlet) = 0;
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
        if ( i > 1)
            GPInfo = EvaluateConstitutiveLaw(GPInfo, XStep, Elements, false, RKMethodLaw);
        end
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod,  dt, false, AlphaStabM);

        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
    
%         invCf = C\(f + uDirichlet);
        k(:,i) = C\(K*XStep + f + uDirichlet);
        if ( loadStep == 1)
            k(:,i) = k(:,i) + (1/dt)* (C\fini);
        end
    end
    XNew  = X; 
   for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end     
    X = XNew;
    
    % Compute stress and D
    GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false, RKMethodLaw);

    if ( drift)
        
        [~,~,~,F] = ComputeForceVector(dt*loadStep, Nodes, Elements, GPInfo, CP);
        UnbalancedForces = fini + F + f0 - ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
        UnbalancedForces(nDirichlet) = 0;
        
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod,  dt, false);
        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
        
        if ( CP.HydroMechanical)
            dX =  (C+dt*K)\(UnbalancedForces);
        else
            dX = C\UnbalancedForces;
        end
        X = X + dX;
        GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false, RKMethodLaw);
        
    end
    
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    
%     PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*loadStep, false, ['ThisProblem-', ElementType]);
    if ( DoSomePostProcess ) 
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
    
end



PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps+0.1, false, ['ThisProblem-', ElementType]);



% Compute the mechanical residual, just to have some fun,....
if (nargout > 2)
    [~,~,~,F] = ComputeForceVector(dt*nSteps, Nodes, Elements, GPInfo, CP);
    % mechanical part
	finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
   
    residual = fini + F + f0 - finter;
    residual(nDirichlet) = 0;
    normResidual = norm(residual);
end
