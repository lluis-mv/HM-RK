function [] = Example0()
close all
addpath('../Sources')
% 1. Define the problem




CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);
[NodesQQ, ElementsQQ] = ReadTheMesh('ThisMeshQ.msh');
[NodesTT, ElementsTT] = ReadTheMesh('ThisMeshT.msh');

[GPInfo] = ComputeElementalMatrices(NodesT, ElementsT, CP, 'T6T3');




DDTT = 10.^linspace(-6,6,60);



% Now the same to compute norms...
for elem = [2, 1]
    figure(30+elem); clf;
    figure(90+elem); clf;
    color = 1;
    
    for Stab = [-1, 0, 1]
        if ( elem == 1)
            ElementType = 'T6T3';
            Nodes = NodesT;
            Elements = ElementsT;
            ThisNumber = 6;
        elseif (elem == 2)
            ElementType = 'Q8Q4';
            Nodes = NodesQ;
            Elements = ElementsQ;
            ThisNumber = 6;
        end
        
        [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
        hhee = [];
        for eell = 1:size(GPInfo, 1)
            hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
        end
        he = mean(hhee);
        
        
        i = 1;
        ddtt = DDTT;
        for dt = DDTT
            
            nSteps = 1;
            
            
            [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
            hhee = [];
            for eell = 1:size(GPInfo, 1)
                hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
            end
            he = mean(hhee);
            if ( Stab < 0)
                %[U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt/nSteps, nSteps, ElementType, 1, Stab);
                [U, GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt/nSteps, nSteps, ElementType, 0);
            else
                [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt/nSteps, nSteps, ElementType, 1, Stab);
            end
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, dt, CP, GPInfo,U);
            
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
            
            
            if ( Stab >= 0)
                [A, nDir, nnoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, 1, Stab);
                nNodes = size(Nodes, 1);
                ii = eye(3*nNodes, 3*nNodes);
                
                B = ii + dt*A;
                B = B(nnoDir, nnoDir);
                
                values = eig(full(B), 'nobalance');
                values = abs(values);
                minval(i)= min(values);
                maxval(i) = max(values);
            end
            
            
            figure(30)
            loglog( ddtt(1:i)/he, L2(1:i), 'k*-.', ddtt(1:i)/he, L2U(1:i), 'rv-.',  ddtt(1:i)/he, LInf(1:i), 'g*-.',  ddtt(1:i)/he, LInfU(1:i), 'bv-.')
            hold on
            xlabel('$\Delta T = c_v \Delta t / h_e^2  $', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
            set(ll, 'interpreter', 'latex')
            drawnow
            hold off
            
            i = i+1;
        end
        figure(30+elem)
        thisColor = 'k';
        NAME = 'Implicit';
        if ( color == 2)
            NAME = 'Explicit. Unstabilized';
            thisColor = 'r';
        elseif ( color == 3)
            NAME = 'Explicit. Stabilized';
            thisColor = 'b';
        end
        color = color + 1;
        
        
        loglog( ddtt, L2U, [thisColor, 'v-.'], 'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| \mathbf{u}-\mathbf{u}_h\|_{L2}$. ', NAME] )
        hold on
        loglog( ddtt, L2, [thisColor, '*-.'],   'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| p_w-p_w^h\|_{L2}$. ', NAME] )
        xlabel('$t$ (s)', 'interpreter', 'latex')
        ylabel('Error norm', 'interpreter', 'latex');
        set(gca, 'FontSize', 13)
        drawnow
        set(gcf, 'Position', [371, 663, 570, 650])
        
        
        ll = legend('location', 'northoutside');
        set(ll, 'interpreter', 'latex')
        if ( Stab >= 0)
            figure(90+elem)
            loglog( ddtt, maxval, [thisColor, 'v-.'], 'MarkerIndices', 1:3:length(ddtt))
            hold on
            ylabel('$\| \lambda \|$', 'interpreter', 'latex')
            xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
            set(gca, 'FontSize', 13)
            drawnow
        end
    end
    figure(30+elem)
    drawnow
    
    yy = ylim();
    xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1]
    plot(xx, yy, 'k-.', 'HandleVisibility','off')
    ylim(yy);
    drawnow
    
    print(['ExampleZero-ErrorNorms-', ElementType], '-dpdf')
    hold off;
    
    figure(90+elem)
    yy = ylim();
    if (yy(1) > 0.1)
        yy(1) = 0.1;
    end
    ylim(yy);
    
    ll = legend('max$(|\lambda|)$. Unstabilized', ...
        'max$(|\lambda|)$. Stabilized', 'location', 'best');
    set(ll, 'interpreter', 'latex')
    drawnow
    yy = ylim();
    yy(1) = 0.1;
    
    xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
    plot(xx, yy, 'k-.', 'HandleVisibility','off')
    ylim(yy);
    drawnow
    print(['ExampleZero-AMatrix-', ElementType], '-dpdf')
    
end




return;


% Now the same to compute norms...
for elem = [2, 1]
    
    
    for Stab = [1, 0]
        figure(30+elem+10*Stab); clf;
        color = 1;
        for nSteps = [10, 100, 1000]
            if ( elem == 1)
                ElementType = 'T6T3';
                Nodes = NodesT;
                Elements = ElementsT;
                ThisNumber = 6;
            elseif (elem == 2)
                ElementType = 'Q8Q4';
                Nodes = NodesQ;
                Elements = ElementsQ;
                ThisNumber = 6;
            end
            [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
            hhee = [];
            for eell = 1:size(GPInfo, 1)
                hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
            end
            he = mean(hhee);
            
            
            i = 1;
            ddtt = DDTT;
            for dt = DDTT
                
                [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt/nSteps, nSteps, ElementType, 1, Stab);
                [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, dt, CP, GPInfo,U);
                
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
                
                
                
                figure(30)
                loglog( ddtt(1:i)/he, L2(1:i), 'k*-.', ddtt(1:i)/he, L2U(1:i), 'rv-.',  ddtt(1:i)/he, LInf(1:i), 'g*-.',  ddtt(1:i)/he, LInfU(1:i), 'bv-.')
                hold on
                xlabel('$\Delta T = c_v \Delta t / h_e^2  $', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex');
                set(gca, 'FontSize', 14)
                drawnow
                ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
                set(ll, 'interpreter', 'latex')
                drawnow
                hold off
                
                i = i+1;
            end
            
            thisColor = 'r';
            if ( color ==  2)
                thisColor = 'b';
            elseif ( color ==  3)
                thisColor = 'g';
            end
            color = color + 1;
            figure(30+elem+10*Stab);
            loglog( ddtt, L2U, [thisColor, 'v-.'], 'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| \mathbf{u}-\mathbf{u}_h\|_{L2}$. $n_{step}$ = ', num2str(nSteps)])
            hold on
            loglog( ddtt, L2, [thisColor, '*-.'],   'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| p_w-p_w^h\|_{L2}$. $n_{step}$ = ', num2str(nSteps)]) ;
            xlabel('t (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 13)
            drawnow
            hold on
            if ( Stab == 0)
                set(gcf, 'Position', [371, 663, 570, 650])
                
                ll = legend('location', 'northoutside');
                set(ll, 'interpreter', 'latex')
            end
            drawnow
        end
        print(['ExampleZero-ErrorNorms-2-', ElementType, '-', num2str(Stab)], '-dpdf')
        hold off;
    end
end


return;



close all; 


% Now the same to compute norms...
for elem = [2, 1]
    
    
    for Stab = [1, 0]
        figure(30+elem+10*Stab); clf;
        color = 1;
        nSteps = 10;
        for RK = [1, 5, 8]
            if ( elem == 1)
                ElementType = 'T6T3';
                Nodes = NodesT;
                Elements = ElementsT;
                ThisNumber = 6;
            elseif (elem == 2)
                ElementType = 'Q8Q4';
                Nodes = NodesQ;
                Elements = ElementsQ;
                ThisNumber = 6;
            end
            [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
            hhee = [];
            for eell = 1:size(GPInfo, 1)
                hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
            end
            he = mean(hhee);
            
            
            i = 1;
            ddtt = DDTT;
            for dt = DDTT
                
                [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt/nSteps, nSteps, ElementType, RK, Stab);
                [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, dt, CP, GPInfo,U);
                
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
                
                
                
                figure(30)
                loglog( ddtt(1:i)/he, L2(1:i), 'k*-.', ddtt(1:i)/he, L2U(1:i), 'rv-.',  ddtt(1:i)/he, LInf(1:i), 'g*-.',  ddtt(1:i)/he, LInfU(1:i), 'bv-.')
                hold on
                xlabel('$\Delta T = c_v \Delta t / h_e^2  $', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex');
                set(gca, 'FontSize', 14)
                drawnow
                ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
                set(ll, 'interpreter', 'latex')
                drawnow
                hold off
                
                i = i+1;
            end
            
            thisColor = 'r';
            if ( color ==  2)
                thisColor = 'b';
            elseif ( color ==  3)
                thisColor = 'g';
            end
            color = color + 1;
            figure(30+elem+10*Stab);
            loglog( ddtt, L2U, [thisColor, 'v-.'], 'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| \mathbf{u}-\mathbf{u}_h\|_{L2}$. RK-', num2str(RK)])
            hold on
            loglog( ddtt, L2, [thisColor, '*-.'],   'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$\| p_w-p_w^h\|_{L2}$. RK-', num2str(RK)]) ;
            xlabel('t (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 13)
            drawnow
            hold on
            if ( Stab == 0)
                set(gcf, 'Position', [371, 663, 570, 650])
                
                ll = legend('location', 'northoutside');
                set(ll, 'interpreter', 'latex')
            end
            drawnow
        end
          yy = ylim();
    
        xx = 10*(he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
        plot(xx, yy, 'k-.', 'HandleVisibility','off')
        ylim(yy);
        drawnow
    
        print(['ExampleZero-ErrorNorms-3-', ElementType, '-', num2str(Stab)], '-dpdf')
        hold off;
    end
end









function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 1, t, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, Ca, Ka);

Aa = Ca\(Ka);

% [vectors, values] = eig(full(Aa), 'nobalance');
[vectors, values] = eig(full(Aa), eye(size(Aa)), 'qz');
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
