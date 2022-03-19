% Solver for a linear problem

function [X, GPInfo, normRes, ThisInfo, nZero] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, AlphaStab)


if (nargout >= 4)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

if ( nargin < 7)
    AlphaStab = 0;
end

nNodes = size(Nodes, 1);
%nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPInfo] = InitializeConstitutiveLaw(CP, GPInfo);


[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 0, dt, true, AlphaStab);
[~, ~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

if ( any([GPInfo.VonMises] == true) )
    addpath('../ModifiedCamClay/vonMises/')
elseif ( any([GPInfo.MCC] == true) )
    addpath('../ModifiedCamClay/')
    addpath('../ModifiedCamClay/Hashiguchi/')
end


[f0] = ComputeInternalForces(Elements, GPInfo, X, CP.HydroMechanical);
f0(nDirichlet) = 0;

fin_n = f0;
fext_n = 0*fini;

t = 0;

[~, uDirichlet] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP);
uDirichlet = 0*uDirichlet;


if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ImplicitProblem-', ElementType]);
end
if ( DoSomePostProcess ) 
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end


for loadStep = 1:nSteps
    
    Xn = X;
    iter = 0;
    
    t = t + dt;
    [df, vDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
    fext_n1 = fext_n + dt*df;
    if ( loadStep == 1)
        fext_n1 = fext_n1+fini;
    end
    uDirichlet = uDirichlet + dt*vDirichlet;
    
    
    
    while( true )
        
        % Compute D, Sigma...
        GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, Xn, Elements, true);
        
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 0, dt, true, AlphaStab);
        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);
        
        
        % mechanical part
        fin_n1 = ComputeInternalForces( Elements, GPInfo, Xn, CP.HydroMechanical);
        
        residual = (fext_n1 - fext_n) + (fin_n - fin_n1);
        
        if ( CP.HydroMechanical)
            
            res2 = C*(Xn-X)-dt*K*Xn;
        
            for jj = 1:nNodes
                residual(3*jj) = -res2(3*jj);
            end
        end
        
        
        residual(nDirichlet) = uDirichlet(nDirichlet)-Xn(nDirichlet);
        
        
        normRes = norm(residual);
        
        
%         disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        if ( iter > 10)
            disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        end
        if ( normRes < 1E-12 && iter > 0)
%             disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
            
            break;
        end
        if ( iter == 30)
            X = nan*X; 
            Xn = nan*X;
            normRes = nan;
            nZero = nnz(A);
            return;
        end
        if ( iter > 25)
            hola = 1;
        end
        
        if ( CP.HydroMechanical)
            A = C-dt*K;
        else
            A = C;
        end
        
        
        dX = A\residual;
        if ( any(isnan(dX)))
            hola = 1;
        end
        
        Xn = Xn + dX;
        
        
        iter = iter+1;
        
%         CheckNumericalDerivative( Nodes, Elements, GPInfo, CP, ElementType, dt, A, X, Xn);
    end
    
    X = Xn;
    
    GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
    
    fin_n = fin_n1;
    
    %PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*loadStep, false, ['ImplicitProblem-', ElementType]);
    if ( DoSomePostProcess ) 
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
end


nZero = nnz(A);
