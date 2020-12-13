% Solver for a linear problem

function [X, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType)
                            

nNodes = size(Nodes, 1);
%nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 0, dt, true);
[~, ~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);



[GPInfo] = InitializeConstitutiveLaw(GPInfo);
[f0] = ComputeInternalForces(Elements, GPInfo, X);
f0(nDirichlet) = 0;

fin_n = f0;
fext_n = fini;

t = 0;

[~, uDirichlet] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP);
uDirichlet = 0*uDirichlet;

for loadStep = 1:nSteps
    
    Xn = X;
    iter = 0;
    
    t = t + dt;
    [df, vDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
    fext_n1 = fext_n + dt*df;
    uDirichlet = uDirichlet + dt*vDirichlet;
    
    
    
    while( true )
        
        % Compute D, Sigma...
        GPInfo = EvaluateConstitutiveLaw(GPInfo, Xn, Elements, true);
        
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 0, dt, true, 1);
        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
        
        
        
        
        % mechanical part
        fin_n1 = ComputeInternalForces( Elements, GPInfo, Xn);
        
        
        residual = (fext_n1 - fext_n) + (fin_n - fin_n1);
        %residual = fini + f*(loadStep/nSteps) + f0 - fin_n1;
        
        % hydraulical part
        res2 = C*(Xn-X)-dt*K*Xn;
        
        for jj = 1:nNodes
            residual(3*jj) = -res2(3*jj);
        end
        
        
        residual(nDirichlet) = uDirichlet(nDirichlet)-Xn(nDirichlet);
        
        
        normRes = norm(residual/nNodes);
        
        if ( iter > -10)
            disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        end
        if ( normRes < 1E-12 && iter > 0)
            break;
        end
        if ( iter == 30)
            break;
        end
        if ( iter > 25)
            hola = 1;
        end
        
        A = C-dt*K;
        dX = A\residual;
        if ( any(isnan(dX)))
            hola = 1;
        end
        
        Xn = Xn + dX;
        
        
        iter = iter+1;
        
        %CheckNumericalDerivative( Nodes, Elements, GPInfo, CP, ElementType, trash, dt, A, X, Xn);
    end
    
    X = Xn;
    
    GPInfo = FinalizeConstitutiveLaw(GPInfo);
    
    fin_n = fin_n1;
    
    
end


function CheckNumericalDerivative( Nodes, Elements, GPInfo, CP, ElementType, trash, dt, A, X, Xn)
delta = 1E-6;
for i = 1:length(Xn)
    X2 = Xn;
    X2(i) = X2(i)+delta*1i;
    ThisResidual = ComputeThisResidual( Nodes, Elements, GPInfo, CP, ElementType, trash, dt, X, X2);
    deri = imag(ThisResidual)/delta;
    J(:,i)= deri;
end

hola = 1;

function residual = ComputeThisResidual( Nodes, Elements, GPInfo, CP, ElementType, trash, dt, X, Xn)



nNodes = size(Nodes, 1);
%nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, trash, dt, true);

[~, ~, X, ~, ~, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);



[GPInfo] = InitializeConstitutiveLaw(GPInfo);
[f0] = ComputeInternalForces(Elements, GPInfo, X);
f0(nDirichlet) = 0;

% Compute D, Sigma...
GPInfo = EvaluateConstitutiveLaw(GPInfo, Xn, Elements, true);


% Create again C with the appropriate ElastoPlastic stiffness matrix
[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, trash, dt, true, 1, 2,0);

[C, K,  ~, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);



% mechanical part
finter = ComputeInternalForces( Elements, GPInfo, Xn);


residual = fini + f*(1/1) + f0 - finter;

% hydraulical part
res2 = C*(Xn-X)-dt*K*Xn;

for jj = 1:nNodes
    residual(3*jj) = -res2(3*jj);
end


residual(nDirichlet) = -Xn(nDirichlet);
