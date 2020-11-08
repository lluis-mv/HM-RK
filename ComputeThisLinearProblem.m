% Solver for a linear problem

function [X, GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod)




nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false);


[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);




[C2, K2 ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, 1, 2, 0);

[C2, ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C2, K2);


PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);


A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii+dt*A;


if ( RKMethod)
    [a,b] = GetRungeKutta(RKMethod);
end


invCf = (C2\f);
invCfini = (C2\fini);

k = zeros(  3*nNodes, length(b));
% X = B*X+invCfini;
for i = 1:length(b)
    XStep = X;
    for j = 1:i-1
        XStep = XStep + dt*a(i,j)*k(:,j);
    end
    k(:,i) = A*XStep + invCfini;
end
XNew = X;
for i = 1:length(b)
    XNew = XNew + dt*b(i)*k(:,i);
end
X = XNew;


% for i = 2:nSteps
%     X = B*X + (1/nSteps)* invCf;
% end

for t = 2:nSteps
    for i = 1:length(b)
        XStep = X;
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
        k(:,i) = A*XStep + (1/nSteps)*invCf;
    end
    XNew = X;
    for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end
    X = XNew;
end



GPInfo = EvaluateConstitutiveLaw(GPInfo, X, Elements, false);
GPInfo = FinalizeConstitutiveLaw(GPInfo);
PostProcessResults(Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);




