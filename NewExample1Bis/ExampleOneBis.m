function [] = ExampleOneBis()
figure(30); clf;
figure(50); clf;
figure(900); clf;

addpath('../')
% 1. Define the problem

T = 1;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1E-3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;
t = 1;



[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);


ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.04, 0.035, 0.03, 0.025];
ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025];


figure(50); clf;


ESIZE = [1:6];



for Elem = [2,1]
    
    
    figure(30); clf;
    figure(50); clf;
    figure(900); clf;
    figure(1); clf;
    color = 1;
    for MyNumber = [1, 10, 100]
        
        esizeAxis = ESIZE;
        i = 1;
        for eSize = ESIZE
            
            
            [NodesQ, ElementsQ] = ReadTheMesh(['Mesh',num2str(eSize),'.msh']);
            [NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);
            
            if ( Elem == 1)
                ElementType = 'T6T3';
                Nodes = NodesT;
                Elements = ElementsT;
                ThisNumber = 6;
            elseif (Elem == 2)
                ElementType = 'Q8Q4';
                Nodes = NodesQ;
                Elements = ElementsQ;
                ThisNumber = 6;
            end
            
            figure(380); clf;
            if ( eSize == 1)
                plotNodes = false;
            else
                plotNodes = false;
            end
            PlotMesh(Nodes, Elements, plotNodes);
            axis off
            axis equal
            ylim([0,1])
            print(['ExampleOneBis-FemMesh-' ElementType, '-' ,num2str(eSize)], '-dpdf')
            
            
            % Estimate the element size
            
            [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
            hhee = [];
            for eell = 1:size(GPInfo, 1)
                hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
            end
            he = mean(hhee);
            
            esizeAxis(i)=he;
            
            dt = he^2/(MyNumber*CP.k*CP.M)
            
            nSteps = ceil(t/dt)
%             nSteps = MyNumber
            dt = t/nSteps;
            dt = 1/nSteps;
            
            RKMethod = 8;
            [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, 1);
            
            PlotNodal(Nodes, Elements, 1000*U(3:3:end), true);
            
            
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, U);
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP, t);
            
            
            %figure(2105); clf; PlotNodal(Nodes, Elements, U(3:3:end)); colorbar; hold off;
            %figure(2106); clf; PlotNodal(Nodes, Elements,Xa(3:3:end)); colorbar; hold off;
            %figure(2107); clf; PlotNodal(Nodes, Elements, U(3:3:end)-Xa(3:3:end)); colorbar; hold off;
            %figure(2108); clf; PlotNodal(Nodes, Elements, U(3:3:end)-Xa(3:3:end), true); colorbar; hold off;
            
            
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
        name = ['$L_2 p_w;$ $\beta$ =', num2str(MyNumber)];
        merda = loglog(esizeAxis, L2, [thisColor, '*-.']);
        set(merda, 'DisplayName', name);
        
        hold on
        name = ['$L_2u;$ $\beta$ =', num2str(MyNumber)];
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
        esizeAxis
        color = color+1;
    end
    
    
end



function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% Other analytical solution...
nNodes = size(Nodes,1);


for nod = 1:nNodes
%     z= 1-Nodes(nod,2);
    z = 1-Nodes(nod,2);
    [u, pw] = EvaluateSticklePastor( z, t, 1, CP.M*CP.k, CP.M);
    
    Xa(3*(nod-1)+2) = -u;
    Xa(3*(nod-1)+3) = pw;
end


if (any(isnan(Xa)))
    Xa = nan*Xa;
end

figure(900)
subplot(2,1,1)
plot(Nodes(:,2), Xa(2:3:end), 'b*', Nodes(:,2), Xnum(2:3:end), 'r*')
hold off
subplot(2,1,2)
plot(Nodes(:,2), Xa(3:3:end), 'b*', Nodes(:,2), Xnum(3:3:end), 'r*')
hold off
