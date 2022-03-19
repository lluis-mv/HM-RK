
function [X, GPElements, GPNodes, normRes, ThisInfo, nZero] = ComputeImplicitNonLinearProblemNodal(Nodes, Elements, CP, dt, nSteps, ElementType, AlphaStab)

if (nargout >= 5)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end
if ( nargin < 7)
    AlphaStab = 0;
end
nNodes = size(Nodes, 1);
%nElements = size(Elements, 1);

[GPElements] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPElements] = InitializeConstitutiveLaw(CP, GPElements);
[GPNodes] = ConstructNodalIntegrationPoints(CP, Nodes, Elements, GPElements);



[C, K] = EnsambleNodalMatrices(Nodes, Elements, GPElements, GPNodes, CP, ElementType, 0, dt, true, AlphaStab);

[~, ~, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPElements, C, K);

if ( any([GPElements.VonMises] == true) )
    addpath('../ModifiedCamClay/vonMises/')
elseif ( any([GPElements.MCC] == true) )
    addpath('../ModifiedCamClay/')
end


[f0] = ComputeNodalInternalForces(Nodes, Elements, GPElements, GPNodes, X);
f0(nDirichlet) = 0;

fin_n = f0;
fext_n = 0*fini;

t = 0;

[~, uDirichlet] = ComputeForceVector(0, Nodes, Elements, GPElements, CP);
uDirichlet = 0*uDirichlet;

% % % % PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPElements, 0, true, ['ImplicitProblem-', ElementType]);
if ( DoSomePostProcess )
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPElements, X, CP);
end


for loadStep = 1:nSteps
    
    Xn = X;
    iter = 0;
    
    t = t + dt;
    [df, vDirichlet] = ComputeForceVector(t, Nodes, Elements, GPElements, CP);
    fext_n1 = fext_n + dt*df;
    if ( loadStep == 1)
        fext_n1 = fext_n1+fini;
    end
    uDirichlet = uDirichlet + dt*vDirichlet;
    
    
    
    while( true )
        
        % Compute D, Sigma...
        [GPElements, GPNodes] = EvaluateConstitutiveLawNodal(CP, GPElements, GPNodes, Xn, true);
        
        % Create again C with the appropriate ElastoPlastic stiffness matrix
        [C, K] = EnsambleNodalMatrices(Nodes, Elements, GPElements, GPNodes, CP, ElementType, 0, dt, true, AlphaStab, C, K);
        [C, K,  ~, ~, ~] = ApplyBoundaryConditions(Nodes, Elements, GPElements, C, K);
        
        
        % mechanical part
        
        fin_n1 = ComputeNodalInternalForces(Nodes, Elements, GPElements, GPNodes, Xn);
        
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
%             disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
        end
        if ( normRes < 1E-12 && iter > 0)
            break;
        end
        if ( iter == 30)
            %             X = nan*X;
            %             Xn = nan*X;
            %             normRes = nan;
            %             return;
            break
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
        
    end
    
    X = Xn;
    
    GPNodes = FinalizeConstitutiveLaw(CP, GPNodes);
    GPElements = ComputeConstrainedModulus(CP, Nodes, Elements, GPElements, GPNodes);
    fin_n = fin_n1;
    
    %PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*loadStep, false, ['ImplicitProblem-', ElementType]);
    if ( DoSomePostProcess )
        ThisInfo = DoThisPostProcess( loadStep*dt, Nodes, Elements, GPElements, X, CP, ThisInfo, fin_n1);
    end
end

nZero = nnz(A);


function [Cn, Kn] = EnsambleNodalMatrices(Nodes, Elements, GPElements, GPNodes, CP, ElementType, RKMethod, dt, implicit, AlphaStabM, Cn, Kn)

nDofs = 3;
nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


Cn = sparse(nDofs*nNodes, nDofs*nNodes);
Kn = sparse(nDofs*nNodes, nDofs*nNodes);

% if ( nargin == 10)
%     Cn = sparse(nDofs*nNodes, nDofs*nNodes);
%     Kn = sparse(nDofs*nNodes, nDofs*nNodes);
% else
%     Cn(1:3:end, 1:3:end) = 0;
%     Cn(1:3:end, 2:3:end) = 0;
%     Cn(2:3:end, 1:3:end) = 0;
%     Cn(2:3:end, 2:3:end) = 0;
%     
%     for nod = 1:nNodes
%         CPatch = GPNodes(nod).NeigNodes;
%         
%         dofsU = [];
%         dofswP = [];
%         for mm = 1:length(CPatch)
%             dofsU = [dofsU, (CPatch(mm)-1)*nDofs+[1, 2]];
%             dofswP = [dofswP, (CPatch(mm)-1)*nDofs+3];
%         end
%         
%         % Adding effective stresses Kuu
%         Cn(dofsU,dofsU) = Cn(dofsU,dofsU) + GPNodes(nod).B'*GPNodes(nod).D*GPNodes(nod).B*GPNodes(nod).Weight;
%     end
%     return;
% end




perme = CP.k;

mIdentity = [1,1,0]';

AlphaStab = 0;

for nod = 1:nNodes
    CPatch = GPNodes(nod).NeigNodes;
    
    dofsU = [];
    dofswP = [];
    for mm = 1:length(CPatch)
        dofsU = [dofsU, (CPatch(mm)-1)*nDofs+[1, 2]];
        dofswP = [dofswP, (CPatch(mm)-1)*nDofs+3];
    end
    
    % Adding effective stresses Kuu
    Cn(dofsU,dofsU) = Cn(dofsU,dofsU) + GPNodes(nod).B'*GPNodes(nod).D*GPNodes(nod).B*GPNodes(nod).Weight;
    
    
    
    % Adding internal forces due to water phase Kuwp
    for el = GPNodes(nod).NeigElement'
        weight = GPElements(el).Weight;
        N = zeros(1, length(dofswP));
        Cel = Elements(el,:);
        for ind = 1:3
            index = find(dofswP == 3*(Cel(ind)-1)+3);
            N(index) = 7/108;
        end
        index = find(dofswP == 3*(nod-1)+3);
        N(index) = 11/54;
        
        Cn(dofsU, dofswP) = Cn(dofsU, dofswP) - GPNodes(nod).B'*mIdentity* N*weight;
        Cn(dofswP, dofsU) = Cn(dofswP, dofsU) + N'*mIdentity'*GPNodes(nod).B*weight;
    end
    
    Kn(dofswP, dofswP) = Kn(dofswP, dofswP) - GPNodes(nod).dN_dX'*perme*GPNodes(nod).dN_dX*GPNodes(nod).Weight;
    
    
end

ngp = 1;
for el = 1:nElements
    ConstModulus=  GPElements(el,ngp).ConstrainedModulus;
    he = sqrt( sum([GPElements(el,:).Weight]));
    AlphaStab = 2/ConstModulus - dt*perme/he^2/1200;
    AlphaStab = max(0.0, AlphaStab);
    if ( length(AlphaStabM) == 1)
        AlphaStab = -AlphaStab*AlphaStabM;
    elseif ( length(AlphaStabM) == 2)
        AlphaStab = AlphaStabM(1)/ConstModulus - dt*perme/he^2*AlphaStabM(2);
        AlphaStab = -max(0.0, AlphaStab);
    end
    
    dofswP = GPElements(el).dofsWP;
    
    Cn(dofswP, dofswP) = Cn(dofswP, dofswP) - GPElements(el).Ms*AlphaStab * GPElements(el).Weight;
end




function [f] = ComputeNodalInternalForces(Nodes, Elements, GPElements, GPNodes, X )
nDofs = 3;

nNodes = size(Nodes, 1);


sign = -1;
Idev = eye(6);

f = zeros(size(X));
mIdentity = [1;1;0];


for nod = 1:nNodes
    
    CPatch = GPNodes(nod).NeigNodes;
    
    dofsU = [];
    dofswP = [];
    for mm = 1:length(CPatch)
        dofsU = [dofsU, (CPatch(mm)-1)*nDofs+[1, 2]];
        dofswP = [dofswP, (CPatch(mm)-1)*nDofs+3];
    end
    
    ThisStress = Idev*GPNodes(nod).StressNew;
    f(dofsU) = f(dofsU) + GPNodes(nod).B'*(  ThisStress([1,2,4]) ) * GPNodes(nod).Weight;
    
    % Adding internal forces due to water phase Kuwp
    for el = GPNodes(nod).NeigElement'
        weight = GPElements(el).Weight;
        N = zeros(1, length(dofswP));
        wPn = zeros(length(dofswP), 1);
        
        Cel = Elements(el,:);
        for ind = 1:3
            index = find(dofswP == 3*(Cel(ind)-1)+3);
            N(index) = 7/108;
            wPn(index) = X( dofswP(index));
        end
        index = find(dofswP == 3*(nod-1)+3);
        N(index) = 11/54;
        wPn(index) = X( dofswP(index));
        wP = N*wPn;
        
        f(dofsU) = f(dofsU) - GPNodes(nod).B'*mIdentity* wP*weight;
    end
    
    
    
    
end






function [GPElements, GPNodes] = EvaluateConstitutiveLawNodal(CP, GPElements, GPNodes,  U, consistent, RKMethod)

nDofs = 3;
if ( nargin == 5)
    RKMethod = 0;
end

nNodes = size(GPNodes, 1);

for nod = 1:nNodes
    
    
    CPatch = GPNodes(nod).NeigNodes;
    
    dofsU = [];
    for mm = 1:length(CPatch)
        dofsU = [dofsU, (CPatch(mm)-1)*nDofs+[1, 2]];
    end
    Uel = U(dofsU);
    
    GPNodes(nod).StrainNew([1,2,4]) = GPNodes(nod).B*Uel;
    
    GPNodes(nod) = EvaluateLaw( CP, GPNodes(nod), consistent, RKMethod);
    
end



function GP = EvaluateLaw(CP, GP, consistent, RKMethod)

if (GP.MCC)
    
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    X = [GP.StressPrev; GP.HistoryPrev];
    
    if (consistent)
        if ( CP.Elastic)
            [Xnew, D, ~] = ExplicitCamClayE(X, DeltaStrain, -1);
        else
            if (CP.MCC == 20)
                [Xnew, D] = ImplicitCASM(X, DeltaStrain);
            else
            [Xnew, D, ~] = ExplicitCamClay2(X, DeltaStrain, -1);
            
            end
        end
    else
        if ( RKMethod == 1)
            RKMethod = RKMethod+1;
        end
        if ( CP.Elastic)
            [Xnew, ~, D] = ExplicitCamClayE(X, DeltaStrain, RKMethod, false);
        else
            [Xnew, ~, D] = ExplicitCamClay2(X, DeltaStrain, RKMethod, false);
        end
    end
    
    
    GP.StressNew  = Xnew(1:6);
    GP.HistoryNew = Xnew(7);
    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
elseif (GP.VonMises)
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    X = [GP.StressPrev];
    
    [Xnew, Dconsist, D] = ExplicitVonMises(X, DeltaStrain, -1);
    if (consistent)
        D = Dconsist;
    end
    
    GP.StressNew  = Xnew(1:6);
    
    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
    
else
    GP.StressNew = GP.StressPrev + GP.D6*(GP.StrainNew - GP.StrainPrev);
end




function [GPElements] = ComputeConstrainedModulus(CP, Nodes, Elements, GPElements, GPNodes)



if ( GPElements(1,1).MCC )
    [kappa, lambda, M, nu] = GetConstitutiveParameters();
    
    for el = 1:size(GPElements,1)
        Celem = Elements(el,:);
        for i = 1:length(Celem)
            pNode(i) = -mean(mean(GPNodes(Celem(i)).StressNew(1:3)));
        end
        p = mean(pNode);
        for gp = 1:size(GPElements,2)
            if ( CP.Elastic)
                K = p/kappa;
            else
                K = p/lambda;
            end
            
            GPElements(el, gp).ConstrainedModulus =  3*K*(1-nu)/(1+nu);
        end
    end
end