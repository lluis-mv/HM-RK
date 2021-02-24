function [] = ExampleOneBis()
figure(30); clf;
figure(50); clf;
figure(900); clf;

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





ESIZE = [0.2, 0.1, 0.075, 0.05];
ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.04, 0.035, 0.03, 0.025];
ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025];


figure(50); clf;




for Elem = [1,2,3]
    
    
    figure(30); clf;
    figure(50); clf;
    figure(900); clf;
    figure(1); clf;
    color = 1;
    for MyNumber = [1, 10, 100, 1000]
        
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
            
            dx = 0.4; dy = 1;
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
            
            
            % First part. compute the eigenvalues
            
            
            % Estimate the element size
            
            mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
            Nodesa = mesha.Nodes';
            Elementsa = mesha.Elements';
            
            figure(1);
            clf;
            triplot(Elementsa, Nodesa(:,1), Nodesa(:,2), 'k');
            drawnow
            axis equal
            axis off
            sizeString = num2str(eSize);
            index = find(sizeString == '.');
            sizeString(index) ='_'; 
            print(['ExampleOneBis-FemMesh-',sizeString, '.pdf'], '-dpdf')
            
            [GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
            he = mean(sqrt([GPInfo(:,:).Weight]));
            esizeAxis(i)=he;
            
            
            
            
            
            dt = he^2/(MyNumber*CP.k*CP.M)
            
            nSteps = ceil(t/dt)
            dt = t/nSteps;
            
            RKMethod = 1;
            [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, 1);
            
            
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, U);
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
            
            
            
            figure(30)
            
            loglog( esizeAxis(1:i), L2(1:i), 'k*-.', esizeAxis(1:i), L2U(1:i), 'rv-.',  esizeAxis(1:i), LInf(1:i), 'g*-.',  esizeAxis(1:i), LInfU(1:i), 'bv-.')
            hold on
            xlabel('$h_e$ (m)', 'interpreter', 'latex')
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
            
            
            
            SlopeInfp = []; SlopeL2p = []; SlopeInfU = []; SlopeL2U = [];
            for ii = 2:i
                SlopeInfp(ii) = log10(LInf(ii)/LInf(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
                SlopeL2p(ii) = log10(L2(ii)/L2(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
                SlopeInfU(ii) = log10(LInfU(ii)/LInfU(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
                SlopeL2U(ii) = log10(L2U(ii)/L2U(ii-1)) / log10(esizeAxis(ii)/esizeAxis(ii-1));
            end
            SlopeInfp
            SlopeL2p
            SlopeInfU
            SlopeL2U
            
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
        
        figure(50)
        name = ['$L_2 p_w;$ $\alpha$ =', num2str(MyNumber)];
        merda = loglog(esizeAxis, L2, [thisColor, '*-.']);
        set(merda, 'DisplayName', name);
        
        hold on
        name = ['$L_2u;$ $\alpha$ =', num2str(MyNumber)];
        merda = loglog(esizeAxis, L2U, [thisColor, 'v-.']);
        set(merda, 'DisplayName', name);
        
        xlabel('$h_e$ (m)', 'interpreter', 'latex');
        ylabel('Error norm', 'interpreter', 'latex');
        set(gca, 'FontSize', 14)
        drawnow
        
        %     ylim(yy);
        ll = legend( 'location', 'best');
        set(ll, 'interpreter', 'latex')
        drawnow
        %         xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
        %         xticks(ticks);
        print(['ExampleOneBis-ErrorNorms-', ElementType], '-dpdf')
        hold on
        
        color = color+1;
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

