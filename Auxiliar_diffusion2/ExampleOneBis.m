function [] = ExampleOneBis()

close all
addpath('../')
% 1. Define the problem

T = 1E-2;

CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;

CP.HydroMechanical = true;

NSteps = 10.^linspace(0, 6, 10);
NSteps = 10.^2;
NSteps = floor(NSteps); NSteps = sort(NSteps);
NStepsRef = 1;


ddtt = t./NSteps;

ticks = [min(ddtt), max(ddtt)];
ticks = floor(log10(ticks));
ticks = 10.^(ticks(1):2:ticks(end));
ticks = [ticks, min(ddtt), max(ddtt)];
ticks = 10.^unique( log10(ticks));
ticks = unique(ticks);


ESIZE = [0.2, 0.1, 0.075, 0.05, 0.025, 0.02];
ESIZE = [0.2, 0.15, 0.1, 0.075, 0.05, 0.04, 0.03, 0.025, 0.02];
figure(50); clf;

for Elem = 1:3
    
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
        
        dx = 3*eSize; dy = 1;
        model = createpde(1);
        
        R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
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
        clf;
        %     triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
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
        
        
        nSteps = NSteps(1);
        
        
        dt = he^2/(1*ThisNumber*CP.k*CP.M)
        
%         dt = he^2/(6000*CP.k*CP.M)
        nSteps = ceil(t/dt)
        dt = t/nSteps;
        
        [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, 1, Stab);
        
        [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Nodes, Elements, ElementType, t, GPInfo, CP);
        
        
        
        figure(30)
        
        loglog( esizeAxis(1:i), L2(1:i), 'k*-.', esizeAxis(1:i), L2U(1:i), 'rv-.',  esizeAxis(1:i), LInf(1:i), 'g*-.',  esizeAxis(1:i), LInfU(1:i), 'bv-.')
        hold on
        xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
        ylabel('Error norm', 'interpreter', 'latex');
        set(gca, 'FontSize', 14)
        drawnow
        yy = ylim();
        xx = (he)^2/(6000*CP.k*CP.M*ThisNumber)*[1,1];
        %         plot(xx, yy, 'k-.')
        if ( yy(2) > 1E20)
            yy(2) = 1E20;
        end
        ylim(yy);
        ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
        set(ll, 'interpreter', 'latex')
        drawnow
        hold off
        
        
        
        Slope1 = []; Slope2 = []; Slope3 = []; Slope4 = [];
        for ii = 2:i
            Slope1(ii) = log10(LInf(ii)/LInf(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
            Slope2(ii) = log10(L2(ii)/L2(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
            Slope3(ii) = log10(LInfU(ii)/LInfU(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
            Slope4(ii) = log10(L2U(ii)/L2U(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
        end
        Slope1
        Slope2
        Slope3
        Slope4
        
        i = i+1;
    end
    
    figure(50)
    
    loglog( esizeAxis, LInf, 'g*-.',  esizeAxis, LInfU, 'bv-.', esizeAxis, L2, 'k*-.', esizeAxis, L2U, 'rv-.' )
    hold on
    xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
    ylabel('Error norm', 'interpreter', 'latex');
    set(gca, 'FontSize', 14)
    drawnow
    yy = ylim();
    xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
    %         plot(xx, yy, 'k-.')
    if ( yy(2) > 1E20)
        yy(2) = 1E20;
    end
    %     ylim(yy);
    ll = legend( '$L_\infty p_w$', '$L_\infty u$', '$L_2 p_w$', '$L_2 u$','location', 'best');
    set(ll, 'interpreter', 'latex')
    drawnow
    %         xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
    %         xticks(ticks);
    print(['ExampleOneBis-ErrorNorms-', ElementType, '-', num2str(Stab)], '-dpdf')
    hold on
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
    thisNorm = norm(Xa(ind)-X(ind));
    if (isnan(thisNorm))
        LInfU = nan;
        break;
    end
    LInfU = max(LInfU, thisNorm);
end



alfa = 2/3; beta = 1/6;
N1 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 1/6;
N2 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 2/3;
N3 = [ 1 - alfa - beta, alfa,  beta];

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

for el = 1:nElements
    Cel = Elements(el,:);
    dofsU = GPInfo(el,1).dofsU;
    dofsWP = GPInfo(el,1).dofsWPreal;
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:length(w)
        [Nu, Np] = GetShapeFunctions( al(gp), be(gp), length(dofsU), length(dofsWP));
        
        L2U = L2U + wA*w(gp)* norm( Nu*(X(dofsU)-Xa(dofsU)));
        L2 = L2 + wA*w(gp)*abs( Np * ( X(dofsWP)-Xa(dofsWP)));
    end
    
end

L2 = L2/sum([GPInfo.Weight]);
L2U = L2U/sum([GPInfo.Weight]);


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
M = CP.M;
k = CP.k;
for nod = 1:nNodes
    xx = 1-Nodes(nod,2);
    TT = M * t*k;
    pw = 0;
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw = pw + 2/aux * sin( aux * xx) * exp( - aux^2 * TT);
    end
    Xa(3*(nod-1)+3) = pw;
end




for nod = 1:nNodes
    z= 1-Nodes(nod,2);
    
    uu = z-1;
    
    for m = 0:100
        
        term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
        uu = uu+term;
    end
    Xa(3*(nod-1)+2) = uu/M;
end


if (any(isnan(Xa)))
    Xa = nan*Xa;
end


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