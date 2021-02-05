function [] = ExampleThreeElastic()



addpath('../')
% 1. Define the problem



CP.HydroMechanical = true;
CP.E = 100;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-3;

CP.Elastic = true;



ESIZE = 0.35;
RKReference = 8;
RKMethods = [8,1:7];



for Elem = [1,2,3]
    
    esizeAxis = ESIZE;
    i = 1;
    for eSize = ESIZE
        
        if (Elem == 1)
            ElementType = 'T3T3';
        elseif (Elem == 2)
            ElementType = 'T6T3';
        else
            ElementType = 'T6T6';
        end
        
        model = createpde(1);
        
        
        R1 = [3,5, 0, 1, 3, 3, 0, 0, 0, 0, -3, -3]';

        g = decsg(R1);
        geometryFromEdges(model, g);
        
        if ( Elem == 1)
            mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        else
            mesh = generateMesh(model, 'Hmax', eSize);
        end
        
        Nodes = mesh.Nodes';
        Elements = mesh.Elements';
        
        % 
        figure(1);
        clf;
        if ( Elem == 1)
            triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
        end
        drawnow
        axis equal
        axis off
        print('ExampleThree-Elastic-FemMesh', '-dpdf')
        
        % Estimate the element size
        
        mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodesa = mesha.Nodes';
        Elementsa = mesha.Elements';
        [GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
        he = mean(sqrt( mean([GPInfo(:,:).Weight])));
        
        Nadim = 20;
        
        NSteps = [4, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8];
        for j = 1:length(NSteps)
            for RK = RKMethods
                nSteps = NSteps(j);
                dt = 0.15/nSteps;
                
                disp(RK)
                tic;
                [U,GPInfo, rrr, information] = ComputeNLProblem2(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1, false);
                timing = toc;
                
                
                ThisInfo(RK,j).t = [information.t];
                ThisInfo(RK,j).F = [information.F];
                
                if ( RK == RKReference)
                    Nadim = ThisInfo(RK,j).F(end);
                end
                
                RES(RK,j) = rrr;
                ddtt(RK,j) = dt;
                N(RK,j) = ThisInfo(RK,j).F(end);
                Time(RK,j) = timing;
                
                if ( j < 2 )
                    continue;
                end
                
                figure(2105); clf
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    end
                    loglog(ddtt(jj,:), RES(jj,:), [merda, '*-.'])
                    hold on
                    
                end
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Residual', 'interpreter', 'latex')
                ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                set(ll, 'interpreter', 'latex')
                grid
                drawnow
                hold off
                
                figure(2106); clf
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    end
                    loglog(ddtt(jj,:), N(jj,:), [merda, '*-.'])
                    hold on
                end
                grid
                drawnow
                hold off
                ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                set(ll, 'interpreter', 'latex')
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Footing load', 'interpreter', 'latex')
                
                figure(2107); clf
                ind = find(N == 0);
                N(ind) = nan;
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    end
                    loglog(ddtt(jj,:), abs(N(jj,:)-Nadim), [merda, '*-.'])
                    hold on
                end
                
                ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                set(ll, 'interpreter', 'latex')
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid
                drawnow
                hold off
                
                   
                figure(2108); clf
                ind = find(Time == 0);
                Time(ind) = nan;
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    end
                    loglog(Time(jj,:), abs(N(jj,:)-Nadim), [merda, '*-.'])
                    hold on
                end
                
                ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                set(ll, 'interpreter', 'latex')
                xlabel('$t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid
                drawnow
                hold off
                
                
            end
      % Printing....
        figure(2105)
        pause(1)
        print(['ExampleThree-Elastic-Residual-', num2str(Elem)], '-dpdf')
        figure(2106)
        pause(1)
        print(['ExampleThree-Elastic-Bearing-', num2str(Elem)], '-dpdf')
        figure(2107)
        pause(1)
        print(['ExampleThree-Elastic-Error-', num2str(Elem)], '-dpdf')
        figure(2108)
        pause(1)
        print(['ExampleThree-Elastic-TimeError-', num2str(Elem)], '-dpdf')
        end
        
        clear Time;
        clear RES;
        clear ddtt;
        clear N;
        clear Time
     
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




function [A, nDirichlet, nNoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, RKMethod, AlphaStabM)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);


[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);


[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

A = C\(K);

nNoDir = [];
for i = 1:3*nNodes
    if ( any(i == nDirichlet))
        % do nothing
    else
        nNoDir = [nNoDir, i];
    end
end
