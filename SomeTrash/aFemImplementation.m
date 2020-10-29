function [] = aFemImplementation()

figure(100); hold off

hElem =0.35;

figure(101); hold off
figure(102); hold off


[RefDisp, RefForce] = FEM(hElem, 100, 28, true);
figure(100)


for consistent = [false, true]
    for nSteps = [15, 20, 25, 40,100]
        figure(100)
        hold off
        plot(-RefDisp, -RefForce, 'k', 'linewidth', 2)
        hold on
        
        figure(101); hold off
        figure(102); hold off
        
        for method = [ (1:8)+20, -1]
            
            [DispI, ForceI, Residuals] = FEM(hElem, nSteps, method, consistent);
            
            SPEC = '*-.';
            if method == -1
                SPEC = 'sk-.';
            end
            
            figure(100)
            plot(-DispI, -ForceI, SPEC)
            hold on
            drawnow
            
            figure(101)
            semilogy(0:(length(Residuals)-1), Residuals, SPEC)
            hold on
            drawnow
            
            figure(102)
            loglog(Residuals(1:end-1), Residuals(2:end), SPEC)
            hold on
            drawnow
            
        end
        figure(100)
        print(['BVP-', num2str(nSteps), '-', num2str(consistent)], '-dpdf')
        figure(101)
        print(['BVP-A-', num2str(nSteps), '-', num2str(consistent)], '-dpdf')
        figure(102)
        print(['BVP-B-', num2str(nSteps), '-', num2str(consistent)], '-dpdf')
    end
end
%legend('Reference','RK1', 'RK2', 'RK3', 'RK4', 'RK5', 'RK6', 'RK7', 'RK8', 'Implicit', 'location', 'best')


function [Displacement, ForceDisplacement, AllResiduals] = FEM( eSize, nSteps, method, consist)
figure(1)

if (nargin == 0)
    eSize = 0.55;
end
model = createpde(1);
R1 = [3,4,-3,3,3,-3,1,1,-1,-1]';
dx = 1; dy = 3;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
figure(1)
pdeplot(model)
drawnow


X = mesh.Nodes';
C = mesh.Elements';


thisNode = find(X(:,1) == max(X(:,1)) & X(:,2) == max(X(:,2))  );

nElem = size(C,1);
nNodes = size(X,1);

GPInfo = CreateGPInfo(nElem);

nodesBottom = find(X(:,2) == 0);
nodesTop = find(X(:,2) == max(X(:,2)));


DirichletZero = [2*nodesBottom; 2*nodesBottom-1; 2*nodesTop-1];
DirichletOne = 2*nodesTop;


DirichletZero = [2*nodesBottom; 2*nodesBottom-1; 2*nodesTop-1];
DirichletOne = [2*nodesTop];

U = zeros(nNodes*2,1);


tolerance = 1E-10;

xLoadVector = 0;
yLoadVector = 0;

deltaDisplacement = -0.60/(nSteps-1);
prescribedDisplacement = 0;

fExternal = CreateExternalForce(X, C);

fExternal = 10*fExternal;
Displacement = 0;
ForceDisplacement = -10;

AllResiduals = [];

Lagrangian = true;

if ( Lagrangian )
    


    dofsRestricted = [DirichletZero; DirichletOne];
    nRestricted = length(dofsRestricted);
    nSystemSize = 2*nNodes;
    M = sparse( nSystemSize, nRestricted);
    
    M( dofsRestricted,1:nRestricted) = eye(nRestricted);
    
    Zeros = sparse(zeros(nRestricted,nRestricted));
    
    Multipliers = zeros(nRestricted, 1);
    
    uaux = zeros( 2*nNodes, 1);
    uaux(DirichletOne) = 1;
    Maux = M'*uaux;
    Maux = sparse(Maux);
end
    

for loadStep = 2:nSteps
    
    iter = 0;
    
    
    
    prescribedDisplacement = prescribedDisplacement + deltaDisplacement;
    Displacement(loadStep) = prescribedDisplacement;
    tic
    
    if ( Lagrangian)
        Multipliers = 0*Multipliers;
    end
    
    AllResiduals = [];
    
    while (true)
        
        [K, f, GPInfo] = ConstructySystem(X, C, U, GPInfo, method, consist);
        if ( Lagrangian)
            

            
            K = [K,M; M', Zeros ];
            
            ForceDisplacement(loadStep) = sum( f(DirichletOne));
            
            residual = f-fExternal + M*Multipliers;
            residual = [residual; M'*U - Maux*prescribedDisplacement];
            
        else
            
        
            dofsRestricted = [DirichletZero; DirichletOne];
            nRestricted = length(dofsRestricted);
        
            K(dofsRestricted,:) = 0;
            K(dofsRestricted, dofsRestricted) = eye(nRestricted);
        
            ForceDisplacement(loadStep) = sum( f(DirichletOne));
        
            residual = f-fExternal;
        
        
            residual(DirichletOne) = U(DirichletOne)-prescribedDisplacement;
            residual(DirichletZero)= U(DirichletZero);
        end
        
        normResidual = norm(residual);
        if ( loadStep == 2)
            AllResiduals = [AllResiduals, normResidual];
        end
        disp([ ' Iteration :: ', num2str(iter), ' RESIDUAL :: ', num2str(normResidual)]);
        if ( normResidual < tolerance && iter > 0)
            GPInfo = FinalizeStep(GPInfo);
            if ( consist == false)
                toc;
                return;
            end
            
            break
        end
        
        dU = -K\residual;
        U = U + dU(1:2*nNodes);
        if (Lagrangian)
            Multipliers = Multipliers + dU(2*nNodes+1:end);
        end
        
        iter = iter+1;
        if (iter > 20)
            toc
            return;
        end
    end
    toc
    
    
    %figure(22)
    %plot(-Displacement, -ForceDisplacement, '*-.')
    %drawnow
end

Displacement = Displacement/3;


function [K, f, GPInfo] = ConstructySystem(X, C, U, GPInfo, method, consist)

nElem = size(C,1);
nNodes = size(X,1);
nDim = 2*nNodes;
K = sparse(nDim,nDim);
f = zeros(2*nNodes,1);

% Compute Current total epsilon
for el = 1:nElem
    Ce = C(el,:);
    xE = X(Ce,:);
    
    nSystem = [2*Ce(1)-1,2*Ce(1), 2*Ce(2)-1, 2*Ce(2), 2*Ce(3)-1, 2*Ce(3)];
    Uel = U(nSystem);
    
    [B, Area] = ComputeBandA(xE(1,1), xE(1,2), xE(2,1), xE(2,2), xE(3,1),xE(3,2) );
    
    GPInfo(el).strainNew = B*Uel;
    
    GPInfo(el) = EvaluateConstitutiveLaw( GPInfo(el) , method, consist);
    
    stress = SetAppropiateSize( GPInfo(el).stressNew );
    
    fE = B'*stress*Area;
    
    f(nSystem) = f(nSystem) + fE;
    
    kE = B'*GPInfo(el).D*B*Area;
    
    K(nSystem,nSystem) = K(nSystem, nSystem) + kE;
    
end



function stress = SetAppropiateSize( stress)
stress = stress([1,2,4]);

function GP = EvaluateConstitutiveLaw( GP , method, consist)

epsilon = zeros(6,1);
epsilon([1,2,4]) = GP.strainNew;
GP.strainNew = epsilon;

DeltaStrain = GP.strainNew - GP.strainPrev;

elastic = false;
if (elastic == true)
    
    
    De = zeros(6,6);
    
    E = 1000;
    nu = 0.3;
    
    for i = 1:3
        for j = 1:3
            if ( i==j)
                De(i,j) = (1-nu);
            else
                De(i,j) = nu;
            end
        end
        De(i+3, i+3) = (1-2*nu)/2;
    end
    
    De = E/(1+nu)/(1-2*nu) * De;
    
    GP.stressNew = GP.stressPrev + De*DeltaStrain;
    GP.D = De([1,2,4], [1,2,4]);
else
    X = [GP.stressPrev; GP.historyPrev];
    
    
    %     [X, DeltaStrain] = FromMechanicsToGeotechnics( X, DeltaStrain);
    
    [Xnew, D, D2] = ExplicitCamClay(X, DeltaStrain, method);
    
    if ( consist == false)
        D = D2;
    end
    
    
    %     [Xnew] = FromGeotechnicsToMechanics( Xnew);
    GP.stressNew  = Xnew(1:6);
    GP.historyNew = Xnew(7);
    GP.D = D([1,2,4], [1,2,4]);
end


function GP = CreateGPInfo(nElem)

GP.strainPrev  = zeros(6,1);
GP.strainNew   = zeros(6,1);
GP.stressPrev  = 10*[ones(3,1); zeros(3,1)];
GP.stressNew   = GP.stressPrev;
GP.historyPrev  = 5;
GP.historyNew  = 0;
GP.D           = zeros(6,6);

for i = 2:nElem
    GP(i) = GP(1);
end

function GP = FinalizeStep(GP)

for i = 1:length(GP)
    GP(i).strainPrev  = GP(i).strainNew;
    GP(i).stressPrev  = GP(i).stressNew;
    GP(i).historyPrev = GP(i).historyNew;
end



function fExt = CreateExternalForce(X, C)



nElem = size(C,1);
nNodes = size(X,1);
fExt = zeros(2*nNodes,1);


nodesLeft = find(X(:,1) == min(X(:,1)) );
nodesRight = find(X(:,1) == max(X(:,1)) );


tLeft = [1;0];
for el = 1:nElem
    Ce = C(el,:);
    ind = [];
    for i = 1:length(Ce)
        ii = find(Ce(i) == nodesLeft);
        if ( length(ii) > 0)
            ind = [ind, i];
        end
    end
    if ( length(ind) ~= 2)
        continue;
    end
    
    xE = X(Ce,:);
    xN = xE(ind,:);
    
    long = norm( xN(2,:)-xN(1,:))/2;
    
    index = Ce(ind);
    tE = eye(2)*tLeft;
    for i = 1:2
        iii = index(i);
        iii = [2*iii-1, 2*iii];
        fExt(iii) = fExt(iii) + tE*long;
    end
end



tRight = -tLeft;
for el = 1:nElem
    Ce = C(el,:);
    ind = [];
    for i = 1:length(Ce)
        ii = find(Ce(i) == nodesRight);
        if ( length(ii) > 0)
            ind = [ind, i];
        end
    end
    if ( length(ind) ~= 2)
        continue;
    end
    
    xE = X(Ce,:);
    xN = xE(ind,:);
    
    long = norm( xN(2,:)-xN(1,:))/2;
    
    index = Ce(ind);
    tE = eye(2)*tRight;
    for i = 1:2
        iii = index(i);
        iii = [2*iii-1, 2*iii];
        fExt(iii) = fExt(iii) + tE*long;
    end
end

