% Solver for a linear problem

function [X, GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

if (nargin ==7)
    AlphaStabM = 1;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);


[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);


PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);


A = C\(K);
ii = eye(3*nNodes, 3*nNodes);


if ( RKMethod)
    [a,b] = GetRungeKutta(RKMethod);
end

C2 = C;
invCf = (C\f);
invCfini = (C\fini);

k = zeros(  3*nNodes, length(b));

for t = 1:nSteps
    for i = 1:length(b)
        XStep = X;
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
         k(:,i) = A*XStep + (1/dt)*invCf;
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




