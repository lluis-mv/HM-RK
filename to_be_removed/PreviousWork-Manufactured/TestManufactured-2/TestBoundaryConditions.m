function [] = ExampleOneBis()
figure(30); clf;
figure(50); clf;
figure(900); clf;

addpath('../')
% 1. Define the problem

T = 1;


CP.HydroMechanical = true;
CP.E = 100;
CP.nu = 0.3;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;
t = 1;



% [NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
% [NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);


ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.04, 0.035, 0.03, 0.025];
ESIZE = [0.2, 0.15, 0.1, 0.075, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025];


figure(50); clf;


ESIZE = [2:5];

figure(50); clf;
color = 1;
for Elem = [1,2]
    
    
    figure(30); clf;
    
    figure(900); clf;
    figure(1); clf;
    
    for MyNumber = [10]
        esizeAxis = ESIZE;
        i = 1;
        for eSize = ESIZE
            
            
            [NodesQ, ElementsQ] = ReadTheMesh(['Mesh',num2str(eSize),'.msh']);
            [NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);
            
            if (Elem == 1)
                ElementType = 'T3T3';
                Nodes = NodesT;
                Elements = ElementsT;
                ThisNumber = 6;
            elseif (Elem == 2)
                ElementType = 'Q4Q4';
                Nodes = NodesQ;
                Elements = ElementsQ;
                ThisNumber = 6;
            end
            [Nodes, Elements] = SimplifyOrder(Nodes,Elements)
            
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
            
            dt = he^2/(MyNumber*CP.k*CP.M);
            
            nSteps = ceil(t/dt);
            %             nSteps = MyNumber
            nSteps = 1;
            dt = t/nSteps;
            
            RKMethod = 1;
            [U,GPInfo] = ComputeLinearProblemDIFFRobin(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, 0);
            
            %PlotNodal(Nodes, Elements, 1000*U(3:3:end), true);
            
            
            
            
            figure(2105); clf; PlotNodal(Nodes, Elements, U(3:3:end)); colorbar; hold off;
            
            caxis([min(U), max(U)])
            axis equal
            xlim([0,0.4])
            ylim([0,1])
            axis off
            pbaspect([1 1 10])
            
            i = i+1;
        end
        
        
    end
    if ( color == 1)
        thisColor = 'k';
    elseif (color == 2)
        thisColor = 'r';
    elseif ( color == 3)
        thisColor = 'g';
    elseif ( color == 4)
        thisColor = 'b';
    end
    
    figure(50)
    name = [ElementType(3:4)];
    merda = loglog(esizeAxis, L2, [thisColor, '*-.']);
    set(merda, 'DisplayName', name);
    
    hold on
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



function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% Other analytical solution...
nNodes = size(Nodes,1);


for nod = 1:nNodes
    %     z= 1-Nodes(nod,2);
    x = Nodes(nod,1);
    y = Nodes(nod,2);
    %     [u, pw] = EvaluateSticklePastor( z, t, 1, CP.M*CP.k, CP.M);
    [u, v, pw] = solution(x, y);
    Xa(3*(nod-1)+1) = u;
    Xa(3*(nod-1)+2) = v;
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
