function [] = ExampleThreeBis()

figure(30); clf;
figure(50); clf;
figure(900); clf;

addpath('../Sources')
% 1. Define the problem

T = 1E-1;

CP.HydroMechanical = false;
CP.E = 1;
CP.nu = 0.49;





ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.04, 0.035, 0.03];
ESIZE = [0.2];


figure(50); clf;
RKMethod = 1;
Elem = 2;

% for Elem = [1, 2, 3]
for Elem = [1,2,3]
    
    esizeAxis = ESIZE;
    i = 1;
    for eSize = ESIZE
        
        if (Elem == 1)
            ElementType = 'T3T3';
            ThisNumber = 200;
        elseif (Elem == 2)
            ElementType = 'T6T3';
            ThisNumber = 6;
        else
            ElementType = 'T6T6';
            ThisNumber = 2000;
        end
        
        dx = 0.3; dy = 1;
        model = createpde(1);
        
        
        R1 = [3,6, 0, 1, 1.35, 4, 4, 0, 0, 0, 0, 0, -4, -4]';
        g = decsg(R1);
        geometryFromEdges(model, g);
        
        if ( Elem == 1)
            mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        else
            mesh = generateMesh(model, 'Hmax', eSize);
        end
        
        Nodes = mesh.Nodes';
        Elements = mesh.Elements';
        
        Stab = 1;
        % First part. compute the eigenvalues
        figure(1);
        pause(1)
        clf;
        if ( Elem == 1)
            triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
        end
        drawnow
        axis equal
        axis off
        print('ExampleOneBis-FemMesh', '-dpdf')
        
        % Estimate the element size
        
        mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodesa = mesha.Nodes';
        Elementsa = mesha.Elements';
        [GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
        he = mean(sqrt( mean([GPInfo(:,:).Weight])));
        esizeAxis(i)=he;
        
        
        
        
        
        nSteps = 10;
        dt = 0.1/nSteps;
        
        
        
        RK = 3;
        
        %[U,GPInfo, information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt/5, 5*nSteps, ElementType);
        [U,GPInfo, trash, information] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, 2, 1, true);
        ThisInfo(Elem).t = [information.t];
        ThisInfo(Elem).F = [information.F];
        
        
        
        figure(233); clf;
 
        plot([ThisInfo(1).t], [ThisInfo(1).F], 'k*-.', 'linewidth', 1.5)
        hold on
        for i = 2:length(ThisInfo)
            figure(233)
            plot([ThisInfo(i).t], [ThisInfo(i).F], '*-.')
            hold on
            drawnow;
        end
        legend('T3T3', 'T6T3', 'T6T6', 'location', 'best')
        
        
        
        
        i = i+1;
        
    end
    
    
end









function [L2, L2U, LInf, LInfU] = ComputeErrorNorms(X, Nodes, Elements, ElementType, t, GPInfo, CP)

[Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, X);

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


indexWP = 3*[1:nNodes];
LInf = max( abs( X(indexWP)-Xa(indexWP)));
LInfU = 0;
for i = 1:nNodes
    ind = 3*(i-1)+[1,2];
    thisNorm = (Xa(ind)-X(ind));
    thisNorm = abs(thisNorm(2));
    if (isnan(thisNorm))
        LInfU = nan;
        break;
    end
    LInfU = max(LInfU, thisNorm);
end




wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;

auxK = [Na2, Na2, wa;
    Na1, Na2, wa;
    Na2, Na1, wa;
    Nb2, Nb2, wb ;
    Nb1, Nb2, wb ;
    Nb2, Nb1, wb ];

al = auxK(:,1)';
be = auxK(:,2)';
w = auxK(:,3)'/sum(auxK(:,3)');


L2 = 0;
L2U = 0;
L2a=0;
L2Ua = 0;
for el = 1:nElements
    Cel = Elements(el,:);
    dofsU = GPInfo(el,1).dofsU;
    dofsWP = GPInfo(el,1).dofsWPreal;
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:length(w)
        [Nu, Np] = GetShapeFunctions( al(gp), be(gp), length(dofsU), length(dofsWP));
        
        L2U = L2U + wA*w(gp)* norm( Nu*(X(dofsU)-Xa(dofsU)))^2;
        L2 = L2 + wA*w(gp)*abs( Np * ( X(dofsWP)-Xa(dofsWP)))^2;
        
        L2Ua = L2Ua + wA*w(gp)*abs( Nu(1,1:2:end)*( X(dofsU(2:2:end))-Xa(dofsU(2:2:end))))^2;
        %         Xpg = Nu(1,1:2:end)*Nodes(Cel,:);
        %         [XX] = ComputeAnalyticalSolution(Xpg, Elements, ElementType, t, CP, GPInfo, X(1:3));
        %         L2U = L2U + wA*w(gp)* norm( Nu*(X(dofsU))-XX(1:2))^2;
        %         L2 = L2 + wA*w(gp)*abs( Np * ( X(dofsWP))-XX(3))^2;
        
    end
    
end

L2 = sqrt(L2/sum([GPInfo.Weight]));
L2U = sqrt(L2U/sum([GPInfo.Weight]));
L2a = sqrt(L2a/sum([GPInfo.Weight]));
L2Ua = sqrt(L2Ua/sum([GPInfo.Weight]));

if (L2 < 1E-15)
    L2 = rand*1E-15;
end
if (L2U < 1E-15)
    L2U = rand*1E-15;
end
if (LInf < 1E-15)
    LInf = rand*1E-15;
end
if (LInfU < 1E-15)
    LInfU = rand*1E-15;
end




function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% Other analytical solution...
nNodes = size(Nodes,1);
E = CP.E;
nu = CP.nu;
for nod = 1:nNodes
    y = Nodes(nod,2);
    
    
    Xa(3*(nod-1)+2) = 0.1*y^2*(y-1)^2*t^2;
    Xa(3*(nod-1)+3) = -(E*t*y*(2*y^2 - 3*y + 1))/(15*(2*nu - 1));
end




figure(900)
subplot(2,1,1)
plot(Nodes(:,2), 0*Xa(2:3:end), 'g*', Nodes(:,2), Xnum(2:3:end), 'r*')
hold off

subplot(2,1,2)
plot(Nodes(:,2), 0*Xa(3:3:end), 'g*', Nodes(:,2), Xnum(3:3:end), 'r*')
hola = 1;
hold off



function [Nu, Np] = GetShapeFunctions( alfa, beta, nU, nP)

ndim = 2;
if (nU == 6)
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
    
    Nu = (zeros(ndim, 3*ndim));
    for i = 1:3
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
elseif (nU == 12)
    Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];
    
    
    Nu = (zeros(ndim, 6*ndim));
    for i = 1:6
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
end


if ( nP == 3)
    Np =  [ 1 - alfa - beta; alfa;  beta]';
elseif ( nP == 6)
    Np =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)]';
end




function [p, y] = CorrectInterpolation(p1, y1)

alfa = linspace(0, 1);
y = [];
p = [];

for ind = 1:2:length(y1)-1
    pp = p1(ind:ind+2);
    yy = y1(ind:ind+2);
    
    N = [(1 - alfa).*(1-2*alfa);
        4*(1-alfa).*alfa;
        alfa.*(2*alfa-1)];
    y = [y; N'*yy];
    p = [p; N'*pp];
end

