function [] = InfluenceStab()
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
% print('ExampleOne-FemMesh', '-dpdf')

Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));


NStepsRef = 10;



if ( true)
    
    
    SStab = 10.^linspace(-2,2,15);
    
    for Elem = [1:3]
        %     for Stab = [1, 0]
        
        
        if ( Elem == 1)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
            ThisNumber = 200;
        elseif (Elem == 2)
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
        
        figure(20);
        hold off;
        
        for RK = [1,3, 8]
            
            firstTime = true;
            i = 1;
            for Stab = SStab
                
                nSteps = NStepsRef;
                dt = t/nSteps;
                [U,GPInfo] = ComputeLinearProblemFast(Nodes, Elements, CP, dt, nSteps, ElementType, RK, -Stab);
                if ( firstTime)
                    [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
                end
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes1, Elements1, GPInfo);
                
                
                
                
% %                 figure(10)
% %                 loglog( SStab(1:i), L2(1:i), 'k*-.', SStab(1:i), L2U(1:i), 'rv-.',  SStab(1:i), LInf(1:i), 'g*-.',  SStab(1:i), LInfU(1:i), 'bv-.')
% %                 hold on
% %                 xlabel('$\beta_s$', 'interpreter', 'latex')
% %                 ylabel('Error norm', 'interpreter', 'latex');
% %                 set(gca, 'FontSize', 14)
% %                 drawnow
% %                 %             yy = ylim();
% %                 %             xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
% %                 %             plot(xx, yy, 'k-.')
% %                 %             if ( yy(2) > 1E20)
% %                 %                 yy(2) = 1E20;
% %                 %             end
% %                 %             ylim(yy);
% %                 ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
% %                 set(ll, 'interpreter', 'latex')
% %                 drawnow
% %                 hold off
                
                i = i+1;
            end
            
            
            
            color = 'r-.';
            if ( RK == 3)
                color = 'b-.';
            elseif (RK == 8)
                color = 'k:';
            end
            
            figure(20)
            loglog(SStab, L2, [color, '*'], 'DisplayName',  ['$L_2 p_w$ RK-', num2str(RK)]);
            hold on
            loglog(SStab, L2U, [color, 'v'], 'DisplayName',  ['$L_2 u$ RK-', num2str(RK)]);
            
            xlabel('Stabilization factor, $\beta_s$', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            ylim([1E-4,1E6])
            xlim([1E-2,1E2])
            legend('location', 'best', 'interpreter', 'latex')
            print(['ExampleOne-Stab-', num2str(Elem)], '-dpdf');
            
        end
        
    end
end



NStepsRef = 5E4;



SStab = [linspace(0, 0.48, 17), linspace(0.5,2,7)];

for Elem = [1:3]
    %     for Stab = [1, 0]
    
    
    if ( Elem == 1)
        ElementType = 'T3T3';
        Nodes = Nodes1;
        Elements = Elements1;
        ThisNumber = 200;
    elseif (Elem == 2)
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
    
    figure(20);
    hold off;
    
    for RK = [1,3, 8]
        
        firstTime = true;
        i = 1;
        for Stab = SStab
            
            nSteps = NStepsRef;
            dt = t/nSteps;
            [U,GPInfo] = ComputeLinearProblemFast(Nodes, Elements, CP, dt, nSteps, ElementType, RK, -Stab);
            if ( firstTime)
                [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
            end
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes1, Elements1, GPInfo);
            
            
            
            
            figure(10)
            semilogy( SStab(1:i), L2(1:i), 'k*-.', SStab(1:i), L2U(1:i), 'rv-.',  SStab(1:i), LInf(1:i), 'g*-.',  SStab(1:i), LInfU(1:i), 'bv-.')
            hold on
            xlabel('$\beta_s$', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            %             yy = ylim();
            %             xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
            %             plot(xx, yy, 'k-.')
            %             if ( yy(2) > 1E20)
            %                 yy(2) = 1E20;
            %             end
            %             ylim(yy);
            ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
            set(ll, 'interpreter', 'latex')
            drawnow
            hold off
            
            i = i+1;
        end
        
        color = 'r-.';
        if ( RK == 3)
            color = 'b-.';
        elseif (RK == 8)
            color = 'k:';
        end
        
        
        
        
        figure(20)
        semilogy(SStab, L2, [color, '*'], 'DisplayName',  ['$L_2 p_w$ RK-', num2str(RK)]);
        hold on
        loglog(SStab, L2U, [color, 'v'], 'DisplayName',  ['$L_2 u$ RK-', num2str(RK)]);
        
        xlabel('Stabilization factor, $\beta_s$', 'interpreter', 'latex')
        ylabel('Error norm', 'interpreter', 'latex');
        set(gca, 'FontSize', 14)
        drawnow
        ylim([1E-14,1])
        legend('location', 'best', 'interpreter', 'latex')
        print(['ExampleOne-Stab2-', num2str(Elem)], '-dpdf');
        
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


