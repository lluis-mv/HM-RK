% Solver for a linear problem

function [X, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps)


nNodes = size(Nodes, 1);
%nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, true);

[~, ~, X, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, C, K);



[GPInfo] = InitializeConstitutiveLaw(GPInfo);


for i = 1:nSteps
    
    Xn = X;
    
    iter = 0;
    while( iter < 30)
        
        % Compute D, Sigma...
        GPInfo = EvaluateConstitutiveLaw(GPInfo, Xn, Elements, true);
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, true, 1, 2,0);

        [C, K,  ~, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);
        
        fini = 2/12*fini;
        
        % mechanical part
        finter = ComputeInternalForces( Elements, GPInfo, Xn);
        
        
        residual = fini + f -finter;
        
        % hydraulical part
        res2 = C*(Xn-X)-dt*K*Xn;
        
        for jj = 1:nNodes
            residual(3*jj) = -res2(3*jj);
        end
        
        
        residual(nDirichlet) = -X(nDirichlet);
        
        
        normRes = norm(residual);
        
        disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        
        if ( normRes < 1E-12 && iter > 0)
            break;
        end
        
        A = C-dt*K;
        dX = A\residual;
        Xn = Xn + dX;
        
        
        iter = iter+1;

    end

    X = Xn;
    %GPInfo = EvaluateConstitutiveLaw(GPInfo, Xn, Elements);
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    
    
end



function [f] = ComputeInternalForces(Elements, GPInfo, X)


nElements = size(Elements, 1);

f = zeros(size(X));
m = [1;1;0];
for el = 1:nElements
    ind = Elements(el,:);
    index = [];
    for ii = 1:length(ind)
        index = [ index, (ind(ii)-1)*3 + [1,2] ];
    end
    
    PW = X( 3*(ind-1)+3);
    
    pw = 1/3*(PW(1)+PW(2)+PW(3));
    
    f(index) = f(index) + GPInfo(el).B'*( GPInfo(el).StressNew([1,2,4]) + m*pw)*GPInfo(el).Weight;
end
