% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

if (nargout == 3)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

if (nargin ==7)
    AlphaStabM = 1;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPInfo] = InitializeConstitutiveLaw(CP, GPInfo);

[C, K ] = AssembleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

fini(nDirichlet) = 0;
A = C\(K);
A(nDirichlet,:)=0;

[~, ~, AllZero] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP);
if ( AllZero)
    invC = 0*C;
else
    invC = inv(C);
end

invCfini = (C\fini);
invCfini(nDirichlet)=0;


if ( RKMethod)
    [a,b, c] = GetRungeKutta(RKMethod);
end
k = zeros(  3*nNodes, length(b));




if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);
end
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
         k(:,i) = A*XStep + invC*(f+uDirichlet);
         if ( loadStep == 1)
             k(:,i) = k(:,i) + (1/dt)*invCfini;
        end
    end
    XNew = X;
    for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end
    X = XNew;
    
    if (any(isnan(X)))
        X = nan*X;
        return;
    end
    if (any(isinf(X)))
        X = nan*X;
        return;
    end
    if ( DoSomePostProcess )
        GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false, dt);
        GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
end



GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false, dt);
GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);
end

