% Solver for a linear problem

function [X, GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

if (nargin ==7)
    AlphaStabM = 1;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[C, K, X, f, fini, nDir] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);


[C2, K2 ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM, 2, 0);
[C2, K2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C2, K2);

PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);

f(nDir) = 0;
fini(nDir) = 0;
A = C\(K);
ii = eye(3*nNodes, 3*nNodes);

A(nDir,nDir)=0;

if ( RKMethod)
    [a,b] = GetRungeKutta(RKMethod);
end

C2 = C;
invCf = (C2\f);
invCfini = (C2\fini);
invCfini(nDir)=0;

k = zeros(  3*nNodes, length(b));

for t = 1:nSteps
    for i = 1:length(b)
        XStep = X;
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
         k(:,i) = A*XStep + (dt/dt)*invCf;
         if ( t == 1)
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
end



GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(GPInfo);
PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);

