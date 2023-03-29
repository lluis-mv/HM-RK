% Solver for a linear problem

function [X, GPInfo, normRes, ThisInfo, nZero] = ComputeImplicitNonLinearProblemContact(Nodes, Elements, CP, dt, nSteps, ElementType, AlphaStab)


Penalty = 10000000;

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

dofsPerNode = 3;
if ( ElementType(1) == 'M')
    dofsPerNode = 4;
end

[~, ~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K, dofsPerNode);

if ( any([GPInfo.VonMises] == true) )
    addpath('../ModifiedCamClay/vonMises/')
elseif ( any([GPInfo.MCC] > 0) )
    addpath('../ModifiedCamClay/')
    addpath('../ModifiedCamClay/VP/')
    addpath('../ModifiedCamClay/Hashiguchi/')
end

    [df, vDirichlet] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP, dofsPerNode);

    dofsPenalty = find(vDirichlet);

    % Remove nDirichelt
    for ii = dofsPenalty'
        index = find(nDirichlet == ii);
        nDirichlet(index) = [];
    end


[f0] = ComputeInternalForces(Elements, GPInfo, X, CP.HydroMechanical);
f0(nDirichlet) = 0;

fin_n = f0;
f0(dofsPenalty) = 0;
fext_n = 0*fini;

t = 0;

[~, uDirichlet] = ComputeForceVector(0, Nodes, Elements, GPInfo, CP, dofsPerNode);
uDirichlet = 0*uDirichlet;


if ( ElementType(1) == 'T')
    PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ImplicitProblem-', ElementType]);
end
if ( DoSomePostProcess ) 
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP, [], [], dofsPerNode);
end

reduce = false;
normRes0 = nan;
proposal = 0*X;

gV = 0.0;

for loadStep = 1:nSteps
    
    Xn = X+ proposal;
    iter = 0;
    
    t = t + dt;
    [df, vDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP, dofsPerNode);

    dofsPenalty = find(vDirichlet);

    % Remove nDirichelt
    for ii = dofsPenalty'
        index = find(nDirichlet == ii);
        nDirichlet(index) = [];
    end


    fext_n1 = fext_n + dt*df;
    if ( loadStep == 1)
        fext_n1 = fext_n1+fini;
    end
    uDirichlet = uDirichlet + dt*vDirichlet;
    
    
    
    while( true )
        
        % Compute D, Sigma...
        GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, Xn, Elements, true, dt);
        
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 0, dt, true, AlphaStab);
        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K, dofsPerNode);
        
        
        % mechanical part
        fin_n1 = ComputeInternalForces( Elements, GPInfo, Xn, CP.HydroMechanical);
        
        residual = (fext_n1 - fext_n) + (fin_n - fin_n1);
        
        if ( CP.HydroMechanical)
            
            res2 = C*(Xn-X)-dt*K*Xn;
            res3 = C*Xn;
            for jj = 1:nNodes
                if ( dofsPerNode == 3)
                    residual(3*jj) = -res2(3*jj);
                elseif (dofsPerNode == 4)
                    residual(4*jj) = -res2(4*jj);
                    residual(4*jj-1) = -res3(4*jj-1);
                end
            end
        end
        
        
        residual(nDirichlet) = uDirichlet(nDirichlet)-Xn(nDirichlet);

        residual(dofsPenalty) = -fin_n1(dofsPenalty) + Penalty * (Xn(dofsPenalty)-gV);

        
        totalForce = sum(Penalty * (Xn(dofsPenalty)-gV));
        residual(length(residual)+1) = totalForce + min(20.0, 10+loadStep);
%         residual(length(residual)+1) = gV-0.0001;


        normRes = norm(residual);
        
        
        disp(['load Step :: ', num2str(loadStep), ' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        if ( iter > 10 || reduce)
            disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
            if ( reduce)
                disp('In the line search, haha')
            end
        end
        if ( normRes < 1E-11 && iter > 0)
            %disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
            
            break;
        end
           if ( normRes < 1E-9 && iter > 11)
            %disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
            
            break;
        end
        if ( iter == 30 && normRes > 1E-8)
            X = nan*X; 
            Xn = nan*X;
            normRes = nan;
            nZero = nnz(A);
            return;
        elseif (iter == 30)
            break;
        end
        if ( iter > 25)
            hola = 1;
        end
        
        if ( CP.HydroMechanical)
            A = C-dt*K;
        else
            A = C;
        end
        

        n =size(A, 1);
        for ii = dofsPenalty'
            A(ii,ii) = A(ii,ii) - Penalty;
            A(ii, n+1) =  - Penalty;
            A(n+1, ii) =  - Penalty;
        end
        A(n+1,n+1) = -length(dofsPenalty)*Penalty;
%         A(n+1,n+1) = 1;

        
        dX = A\residual;
        gV = gV-dX(end);
        dX = dX(1:n);
        if ( any(isnan(dX)))
            hola = 1;
            dX = 0*Xn;
            Xn = X;
            reduce = true;
            iter = 0;
        end
        if ( reduce == true && iter < 10)
            dX = 0.01*dX;
        end
        if (iter > 10 && normRes > normRes0)
            dX = 0.1*(rand()-0.5)*dX;
        end
        normRes0 = normRes;


        Xn = Xn + dX;
        
        
        iter = iter+1;
        
        %B = CheckNumericalDerivative( Nodes, Elements, GPInfo, CP, ElementType, dt, A, X, Xn, AlphaStab, dofsPerNode);
    end
    proposal = Xn-X;
    reduce = false;
    X = Xn;
    
    GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
    
    fin_n = fin_n1;
    
    %PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*loadStep, false, ['ImplicitProblem-', ElementType]);
    if ( DoSomePostProcess ) 
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo, [], dofsPerNode);
    end
end


nZero = nnz(A);
