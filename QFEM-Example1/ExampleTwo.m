function [] = ExampleTwo()
figure(30); clf;
figure(50); clf;
figure(900); clf;

addpath('../Sources')
% 1. Define the problem

T = 1E-2;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


DT = 10.^linspace(-7,2,40);



figure(50); clf;



for Stab = [0,1]
    
    figure(50); clf;
    figure(51); clf;
    color = 1;
    
    
    
    
    for Elem = [1,2,3]
        i = 1;
        
        
        if (Elem == 1)
            ElementType = 'Q4Q4';
        elseif (Elem == 2)
            ElementType = 'Q4Q4';
        else
            ElementType = 'Q8Q4';
        end
        
        dx = 0.4; dy = 1;
        model = createpde(1);
        
        R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
        g = decsg(R1);
        geometryFromEdges(model, g);
        


        [Nodes2, Elements2] = ReadTheMesh('ThisMesh.msh');
        [Nodes, Elements] = SimplifyOrder(Nodes2,Elements2);
        [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'Q4Q4');
        he = mean(sqrt([GPInfo(:,:).Weight]));
        if ( Elem == 3)
            Nodes = Nodes2; 
            Elements = Elements2;
        end
        
        
        
        
        
        
        
        for dt = [DT]
            
            nSteps = 1;
            
            if (Elem ~= 1)
                [U, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, Stab);
            else
                [U, GPInfo] = ComputeImplicitNonLinearProblemNodalQuad(Nodes, Elements, CP, dt, nSteps, ElementType, Stab);
            end
            
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, dt*nSteps, CP, GPInfo, U);
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
            
            
            
            figure(30)
            
            loglog( DT(1:i), L2(1:i), 'k*-.', DT(1:i), L2U(1:i), 'rv-.',  DT(1:i), LInf(1:i), 'g*-.',  DT(1:i), LInfU(1:i), 'bv-.')
            hold on
            xlabel('$dt$ (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            
            ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
            set(ll, 'interpreter', 'latex')
            drawnow
            hold off
            
            DTAxis(i) = dt * CP.M*CP.k/(he^2);
            i = i+1;
        end
        if ( color == 1)
            thisColor = 'r';
        elseif (color == 2)
            thisColor = 'g';
        elseif ( color == 3)
            thisColor = 'b';
        elseif ( color == 4)
            thisColor = 'm';
        end
        ElementType1 = ElementType;
        if ( Elem == 1)
            ElementType1 = ['NS-', ElementType1];
        end
        if (Elem < 3 && Stab > 0)
            ElementType1 = [ElementType1, '. Stab'];
        end
        figure(50)
        semilogx(DTAxis, LInf, [thisColor, 'v-.'], 'DisplayName', ElementType1)
        hold on
        ylim([0, 1])
        %         plot([6,6], [0, 1.2], 'k:', 'HandleVisibility','off')
        ll = legend();
        set(ll, 'interpreter', 'latex', 'location', 'best')
        drawnow;
        hold on;
        color = color+1;
        xlabel('$\Delta t \, c_v / h_e^2$', 'interpreter', 'latex')
        ylabel('$L_\infty p_w$ (kPa)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15);
        
        print(['Qed-LInf', num2str(Stab)], '-dpdf')
        
        figure(51)
        semilogx(DTAxis, L2, [thisColor, 'v-.'], 'DisplayName', ElementType1)
        hold on
        %         plot([6,6], [0, 1.2], 'k:', 'HandleVisibility','off')
        ll = legend();
        set(ll, 'interpreter', 'latex', 'location', 'best')
        drawnow;
        hold on;
        ylim([0, 0.25])
        xlabel('$\Delta t \, c_v / h_e^2$', 'interpreter', 'latex')
        ylabel('$L_2 u$ (m)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15);
        print(['Qed-L2', num2str(Stab)], '-dpdf')
    end
    
    
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

figure(900)
subplot(2,1,1)
plot(Nodes(:,2), Xa(2:3:end), 'b*', Nodes(:,2), Xnum(2:3:end), 'r*')

subplot(2,1,2)
plot(Nodes(:,2), Xa(3:3:end), 'b*', Nodes(:,2), Xnum(3:3:end), 'r*')
hola = 1;

