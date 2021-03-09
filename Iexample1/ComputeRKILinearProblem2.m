% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeRKILinearProblem2(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

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
    
    for iter = 0:1000
        kPrev = k;
        trick = false;
        
        for i = 1:length(b)
            
            prevNumber = 0;
            
            
            XStep = X;
            t = (loadStep-1)*dt + c(i)*dt;
            [f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
            f(nDirichlet) = 0;
            for j = 1:length(b)
                XStep = XStep + dt*a(i,j)*kPrev(:,j);
            end
            
            update =  A*XStep + invC*(f+uDirichlet);
            kNew(:,i) = update;
            thisNorm = norm(kPrev(:,i)-kNew(:,i));
            
            if ( norm(thisNorm) > 10)
                thisNumber = min(10*(1+iter/10)/(norm(update)),  1);
                if ( thisNumber < prevNumber)
                    thisNumber = prevNumber + 0.01;
                    prevNumber = thisNumber;
                else
                    prevNumber = thisNumber;
                end
                update = thisNumber*update;
                trick = true;
            end
            k(:,i) = update;
            
            if ( loadStep == 1)
                k(:,i) = k(:,i) + (1/dt)*invCfini;
            end
            
        end
        
        thisNorm = norm(k-kNew);
        if ( thisNorm < 1E-14 && trick == false)
            break;
        end
        if ( thisNorm < 1E-10 && iter > 100)
            break;
        end
        if ( iter > 900)
            X = X*nan;
            return;
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
PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisProblem-', ElementType]);

