% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeLinearProblemDIFF(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

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

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

RealDofs = unique([unique( [GPInfo.dofsWP]), unique( [GPInfo.dofsU])]);
for i = 1:3*nNodes
    if ( any(i == RealDofs) == false)
        C(i,:) = 0;
        C(:,i) = 0;
        K(i,:) = 0;
        K(:,i) = 0;
        C(i,i) = 1;
        K(i,i) = 1;
    end
end



t = 1;
[f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
f(nDirichlet) = 0;
f = f(3:3:end);
K(nDirichlet, nDirichlet) = eye( length(nDirichlet), length(nDirichlet));
K = K(3:3:end,3:3:end);
U = -K\f;
X(3:3:end) = U;
return;
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
        GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false);
        GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
end



GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);
end

