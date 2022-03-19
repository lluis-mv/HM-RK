
% Solver for a linear problem

function [X, GPInfo, normResidual, ThisInfo] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM, Iterate)

if ( nargin < 9)
    Iterate = true;
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
[GPInfo] = InitializeConstitutiveLaw(CP, GPInfo);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[~,~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

f0 = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
f0(nDirichlet) = 0;

if ( RKMethod)
    [a,b,c] = GetRungeKutta(RKMethod);
end
k = zeros(  3*nNodes, length(b));

if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);
end
if ( DoSomePostProcess )
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end

for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el,gp).DPrev = rand(7,6);
    end
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
            k(:,i) = k(:,i-1);
        end
        kPrev = k(:,i);
        ss = 0;
        while (true)
            GPInfo = EvaluateConstitutiveLawNL(GPInfo, CP,  dt, k, a, b,  i, true);
        
            % Create again C with the appropriate ElastoPlastic stiffness matrix
            [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod,  dt, false, AlphaStabM);
            
            [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
            k(:,i) = C\(K*XStep + f + uDirichlet);
            if ( loadStep == 1)
                k(:,i) = k(:,i) + (1/dt)* (C\fini);
            end
%             norm(k(:,i)-kPrev)
            if (norm(k(:,i)-kPrev) < 1E-12)
                break;
            end
            
            if ( ss > 10)
                disp('breaking the rules')
                break;
            end
            
            if ( Iterate == false && loadStep > 1)
                break;
            end
            kPrev = k(:,i);
            ss = ss+ 1;
        end
            if ( Iterate == false && loadStep > 1)
                kConst = k; kConst(:,i) = kPrev;
                GPInfo = EvaluateConstitutiveLawNL(GPInfo, CP, dt, kConst, a, b, i, false);
            else
                GPInfo = EvaluateConstitutiveLawNL(GPInfo, CP, dt, k, a, b, i, false);
            end
    end
    
    GPInfo = EvaluateConstitutiveLawNL(GPInfo, CP, dt, k, a, b);
    
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


if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps+0.1, false, ['ThisProblem-', ElementType]);
end



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
