
function [X, GPElements, GPNodes, normRes, ThisInfo, nZero] = ComputeImplicitNonLinearProblemNodalQuad(Nodes, Elements, CP, dt, nSteps, ElementType, AlphaStab)

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
GPElements = CalculateSmoothingPathAreas( Nodes, Elements, GPElements);
[GPNodes] = ConstructNodalIntegrationPointsQuad(CP, Nodes, Elements, GPElements);



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

reduce = false;
normRes0 = nan;
proposal = 0*X;
for loadStep = 1:nSteps

    Xn = X+ proposal;
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
        if ( iter > 10 || reduce)
            disp([' :: nonlinear solver, iter :: ', num2str(iter), ' :: residual ', num2str(normRes) ])
            if ( reduce)
                disp('In the line search, haha')
            end
        end
        if ( normRes < 1E-12 && iter > 0)
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


        dX = A\residual;
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

    end
    proposal = Xn-X;
    reduce = false;
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
    Cn(dofsU,dofswP) = Cn(dofsU,dofswP) - GPNodes(nod).Q*GPNodes(nod).Weight;
    Cn(dofswP,dofsU) = Cn(dofswP,dofsU) + GPNodes(nod).Q'*GPNodes(nod).Weight;

    Kn(dofswP, dofswP) = Kn(dofswP, dofswP) - GPNodes(nod).dN_dX'*perme*GPNodes(nod).dN_dX*GPNodes(nod).Weight;


end




ngp = 1;
for el = 1:nElements
    ConstModulus=  GPElements(el,ngp).ConstrainedModulus;
    he = sqrt( sum([GPElements(el,:).Weight]));
    AlphaStab = 1.0/ConstModulus - dt*perme/he^2/120;
    AlphaStab = max(0.0, AlphaStab);
    if ( length(AlphaStabM) == 1)
        AlphaStab = -AlphaStab*AlphaStabM;
    elseif ( length(AlphaStabM) == 2)
        AlphaStab = AlphaStabM(1)/ConstModulus - dt*perme/he^2*AlphaStabM(2);
        AlphaStab = -max(0.0, AlphaStab);
    end


    dofswP = GPElements(el,1).dofsWP;

    Cn(dofswP, dofswP) = Cn(dofswP, dofswP) - GPElements(el,1).Ms*AlphaStab * sum([GPElements(el,:).Weight]);
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
    f(dofsU) = f(dofsU) - GPNodes(nod).Q * X(dofswP) * GPNodes(nod).Weight;
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
