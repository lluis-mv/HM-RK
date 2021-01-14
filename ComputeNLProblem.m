
% Solver for a linear problem

function [X, GPInfo, normResidual, ThisInfo] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

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

% I should get the correct initial D...





for loadStep = 1:nSteps
    for i = 1:length(b)
        
        for sub = [1:5]
            initialize = true;
            if ( sub == 5)
                initialize = false;
            end
            
            XStep = X;
            t = (loadStep-1)*dt + c(i)*dt;
            [f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
            f(nDirichlet) = 0;
            for j = 1:i-1
                XStep = XStep + dt*a(i,j)*k(:,j);
            end
            
            
            % Create again C with the appropriate ElastoPlastic stiffness matrix
            [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod,  dt, false, AlphaStabM);
            
            [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
            
            %         invCf = C\(f + uDirichlet);
            k(:,i) = C\(K*XStep + f + uDirichlet);
            if ( loadStep == 1)
                k(:,i) = k(:,i) + (1/dt)* (C\fini);
            end
            kk(:,sub) = k(:,i);
            if ( sub > 1)
%                 disp(sub)
%                 disp(max( abs(kk(:,sub)-kk(:,sub-1))))
                if ( max( abs(kk(:,sub)-kk(:,sub-1))) == 0)
                    initialize = false;
                end
            end
            
            GPInfo = EvaluateConstitutiveLawNL(GPInfo, X, dt, k, a, b, c, i, initialize);
            if ( i < length(b) && initialize == false)
                GPInfo = EvaluateConstitutiveLawNL(GPInfo, X, dt, k, a, b, c, i+1);
            end
            if ( initialize == false)
                break
            end
        end
        
    end
    
    GPInfo = EvaluateConstitutiveLawNL(GPInfo, X, dt, k, a, b, c);
    
    XNew  = X;
    for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end
    X = XNew;
    
    
    
    %PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*loadStep, false, ['ThisProblem-', ElementType]);
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
    %     figure(50)
    %     triplot(Elements, Nodes(:,1), Nodes(:,2))
    %     axis off; axis equal
    %     hold on
    %     norm(residual)
    %     quiver( Nodes(:,1), Nodes(:,2), residual(1:3:end), residual(2:3:end));
    %     hold off
    %     hola = 1;
end
