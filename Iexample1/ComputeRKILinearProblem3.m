% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeRKILinearProblem3(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

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
    [a,b] = GetRungeKuttaImplicit(RKMethod);
    c = 0*b;
end
k = zeros(  3*nNodes, length(b));




PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);
if ( DoSomePostProcess )
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end

for loadStep = 1:nSteps
    
    k = ComputeTheseK(X,  k, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC);
    
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
PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);




function [k2] = ComputeTheseK(X,  k, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC)


fun = @(x) ComputeKError(x,  X, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC);
options = optimoptions('fsolve', 'FunctionTolerance', 1E-14, 'StepTolerance', 1E-14);
k2 = fsolve(fun, k,options);
k2 = fsolve(fun, k2,options);
% k2 = fsolve(fun, k2,options);
% k2 = fminsearch(fun, k);
hola = 1;



function err = ComputeKError(k,  X, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC)
k2 = ComputeSource(k,  X, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC);
err = k2-k;
norm(err);
% err = norm(err);




function [k2] = ComputeSource(k,  X, a, b, c, dt, Nodes, Elements, GPInfo, CP, nDirichlet, A, invC)

k2 = k; 
for i = 1:length(b)
    
    
    XStep = X;
    t = dt;
    [f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
    f(nDirichlet) = 0;
    for j = 1:length(b)
        XStep = XStep + dt*a(i,j)*k(:,j);
    end
    
    k2(:,i) =  A*XStep + invC*(f+uDirichlet);
    
end



