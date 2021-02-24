function [] = ExampleOne()
close all
addpath('../')
% 1. Define the problem

T = 1E-2;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.075;

model = createpde(1);

dx = 0.4; dy = 1;

R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';


% First part. compute the eigenvalues
figure(1);
clf;
triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
drawnow
axis equal
axis off
print('ExampleOne-FemMesh', '-dpdf')

Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));

NSteps = 10.^linspace(0, 6, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);
NStepsRef = 1;


ddtt = t./NSteps;

ticks = [min(ddtt), max(ddtt)];
ticks = floor(log10(ticks));
ticks = 10.^(ticks(1):2:ticks(end));
ticks = [ticks, min(ddtt), max(ddtt)];
ticks = 10.^unique( log10(ticks));
ticks = unique(ticks);



for Stab = [1, 0]
    firstTime = true;
    for j = 3:-1:1
        
        if ( j == 1)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
        elseif (j == 2)
            ElementType = 'T6T3';
            Nodes = Nodes2;
            Elements = Elements2;
        else
            ElementType = 'T6T6';
            Nodes = Nodes2;
            Elements = Elements2;
        end
        
        
        
        
        
        
        dt = t/NStepsRef;
        [U,GPInfo] = ComputeThisLinearProblemFast(Nodes, Elements, CP, dt,NStepsRef, ElementType, 1, Stab);
        if ( firstTime)
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
            firstTime = false;
            
            index = find(Nodes(:,1) == 0);
            WP = Xa(3*(index-1)+3);
            
            y = Nodes(index,2);
            [ya, index] = sort(y);
            WPa = WP(index);
            figure(987+Stab);
            clf;
            [WPa, ya] = CorrectInterpolation(WPa, ya);
            plot(WPa, ya, 'k', 'linewidth', 2, 'DisplayName', 'Reference')
        end
        
        index = find(Nodes(:,1) == 0);
        WP = U(3*(index-1)+3);
        
        y = Nodes(index,2);
        [y, index] = sort(y);
        WP = WP(index);
        
        if ( all('T6T6' == ElementType))
            [WP, y] = CorrectInterpolation(WP, y);
        end
        
        figure(987+Stab)
        hold on
        if ( all('T3T3' == ElementType))
            plot(WP, y, '-.','DisplayName', ElementType, 'linewidth', 1.5)
        else
            plot(WP, y, 'DisplayName', ElementType, 'linewidth', 1.5)
        end
        ll = legend('location', 'best');
        set(ll, 'interpreter', 'latex')
        xlabel('$p_w$ (kPa)', 'interpreter', 'latex')
        ylabel('$z$ (m)', 'interpreter', 'latex')
        set(gca, 'FontSize', 14)
        print(['ExampleOne-Solution-', num2str(Stab)], '-dpdf')
        
    end
end






if (false)
    
    % First part. Spectral radii
    for j = 1:3
        
        if ( j == 1)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
            ThisNumber = 200;
        elseif (j == 2)
            ElementType = 'T6T3';
            Nodes = Nodes2;
            Elements = Elements2;
            ThisNumber = 6;
        else
            ElementType = 'T6T6';
            Nodes = Nodes2;
            Elements = Elements2;
            ThisNumber = 2000;
        end
        
        
        minval = nan*ddtt;
        maxval = nan*ddtt;
        minval2 = nan*ddtt;
        maxval2 = nan*ddtt;
        i = 1;
        
        for dt = ddtt
            
            [A, nDir, nnoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, 1, 1);
            nNodes = size(Nodes, 1);
            ii = eye(3*nNodes, 3*nNodes);
            
            B = ii + dt*A;
            B = B(nnoDir, nnoDir);
            
            values = eig(full(B), 'nobalance');
            values = abs(values);
            minval(i)= min(values);
            maxval(i) = max(values);
            
            
            
            
            [A, nDir, nnoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, 1, 0);
            nNodes = size(Nodes, 1);
            ii = eye(3*nNodes, 3*nNodes);
            
            B = ii + dt*A;
            B = B(nnoDir, nnoDir);
            
            
            values = eig(full(B), 'nobalance');
            values = abs(values);
            minval2(i)= min(values);
            maxval2(i) = max(values);
            i = i+1;
            
            
            
        end
        figure(j+1);
        clf;
        
        loglog(ddtt, minval2, 'm*-.', ddtt, maxval2, 'c*-.')
        
        drawnow;
        hold on
        loglog(ddtt, minval, 'r*-.', ddtt, maxval, 'b*-.')
        drawnow
        
        xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
        ylabel('$\| \lambda \|$', 'interpreter', 'latex')
        set(gca, 'FontSize', 14)
        drawnow
        yy = ylim();
        xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
        plot(xx, yy, 'k-.')
        ylim(yy);
        ll = legend('min$(|\lambda|)$ Primal', 'max$(|\lambda|)$ Primal', ...
            'min$(|\lambda|)$ Stab', 'max$(|\lambda|)$ Stab', 'location', 'best');
        set(ll, 'interpreter', 'latex')
        xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
        xticks(ticks);
        print(['ExampleOne-Radii-', ElementType], '-dpdf')
    end
    
    
    
    
    
    
    % Now the same to compute norms...
    for Stab = [1, 0]
        for j = 1:3
            
            if ( j == 1)
                ElementType = 'T3T3';
                Nodes = Nodes1;
                Elements = Elements1;
                ThisNumber = 200;
            elseif (j == 2)
                ElementType = 'T6T3';
                Nodes = Nodes2;
                Elements = Elements2;
                ThisNumber = 6;
            else
                ElementType = 'T6T6';
                Nodes = Nodes2;
                Elements = Elements2;
                ThisNumber = 2000;
            end
            
            firstTime = true;
            i = 1;
            for nSteps = NSteps
                
                dt = t/nSteps;
                [U,GPInfo] = ComputeThisLinearProblemFast(Nodes, Elements, CP, dt, nSteps, ElementType, 1, Stab);
                if ( firstTime)
                    [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
                end
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes1, Elements1, GPInfo);
                
                
                
                figure(30+Stab*10+j)
                
                loglog( ddtt(1:i), L2(1:i), 'k*-.', ddtt(1:i), L2U(1:i), 'rv-.',  ddtt(1:i), LInf(1:i), 'g*-.',  ddtt(1:i), LInfU(1:i), 'bv-.')
                hold on
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex');
                set(gca, 'FontSize', 14)
                drawnow
                yy = ylim();
                xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
                plot(xx, yy, 'k-.')
                if ( yy(2) > 1E20)
                    yy(2) = 1E20;
                    yy(1) = 1E-10;
                end
                ylim(yy);
                ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
                set(ll, 'interpreter', 'latex')
                drawnow
                hold off
                
                i = i+1;
            end
            
            figure(30+Stab*10+j)
            clf;
            loglog( ddtt, LInf, 'g*-.',  ddtt, LInfU, 'bv-.', ddtt, L2, 'k*-.', ddtt, L2U, 'rv-.' )
            hold on
            xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            yy = ylim();
            xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
            plot(xx, yy, 'k-.')
            if ( yy(2) > 1E20)
                yy(2) = 1E20;
                yy(1) = 1E-10;
            end
            ylim(yy);
            ll = legend( '$L_\infty p_w$', '$L_\infty u$', '$L_2 p_w$', '$L_2 u$','location', 'best');
            set(ll, 'interpreter', 'latex')
            drawnow
            xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
            xticks(ticks);
            print(['ExampleOne-ErrorNorms-', ElementType, '-', num2str(Stab)], '-dpdf')
            hold off
        end
    end
    
end


if ( true)
    
    
    % Now lets check the RK methods
    
    for j = 1:2
        
        if ( j == 1)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
            ThisNumber = 200;
        elseif (j == 2)
            ElementType = 'T6T3';
            Nodes = Nodes2;
            Elements = Elements2;
            ThisNumber = 6;
        else
            ElementType = 'T6T6';
            Nodes = Nodes2;
            Elements = Elements2;
            ThisNumber = 2000;
        end
        
        Stab = 1;
        
        for RK = [1,3,8]
            firstTime = true;
            i = 1;
            
            for nSteps = NSteps
                
                dt = t/nSteps;
                [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, Stab);
                if ( firstTime)
                    [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
                end
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes1, Elements1, GPInfo);
                
                
                
                i = i+1;
            end
            
            
            color = 'r';
            if ( RK == 3)
                color = 'b';
            elseif (RK == 8)
                color = 'k';
            end
            
            
            figure(20+3*j)
            if (RK == 1)
                clf;
            end
            loglog( ddtt, L2, [color, '*-.'], 'DisplayName',  ['$L_2 p_w$ RK-', num2str(RK)]);
            hold on
            loglog( ddtt, L2U, [color, 'v-.'], 'DisplayName',  ['$L_2 u$ RK-', num2str(RK)]);
            hold on
            xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            
            
        end
        figure(20+3*j)
        
        drawnow
        yy = ylim();
        xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
        plot(xx, yy, 'k-.', 'HandleVisibility','off')
        if ( yy(2) > 1E20)
            yy(2) = 1E20;
            yy(1) = 1E-10;
        end
        ylim(yy);
        drawnow
        xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
        xticks(ticks);
        legend('location', 'best', 'interpreter', 'latex')
        print(['ExampleOne-RK-', ElementType], '-dpdf')
        
        
       
        
    end
end








function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 3, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


Xa = real(Xa);





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
