function [] = ExampleCompressibility()
close all
addpath('../')
% 1. Define the problem




CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);


[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);







DDTT = 10.^linspace(-6,6,60);





% Now the same to compute norms...
for elem = [2, 1]
    
    
    figure(90+elem);
    clf;
    newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    colororder(newcolors)
    
    for Stab = [0]
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
        
        
        
        K = CP.E/3/(1-2*nu);
        for C = [2,10, 100, 1000, 0]*K
            
            CP.Compressibility = 1/C;
            if ( C == 0)
                CP.Compressibility = 0;
            end
            C = C/K;
            i = 1;
            ddtt = DDTT;
            for dt = DDTT
                
                
                
                [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
                hhee = [];
                for eell = 1:size(GPInfo, 1)
                    hhee(eell) =  sqrt( sum([GPInfo(eell,:).Weight]));
                end
                he = mean(hhee);
                
                [A, nDir, nnoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, 1, Stab);
                nNodes = size(Nodes, 1);
                ii = eye(3*nNodes, 3*nNodes);
                
                B = ii + dt*A;
                B = B(nnoDir, nnoDir);
                
                values = eig(full(B), 'nobalance');
                values = abs(values);
                minval(i)= min(values);
                maxval(i) = max(values);
                
                i = i+1;
            end
            
            
            
            figure(90+elem)
            this = num2str(C);
            if ( C == 0)
                this = ' $\infty$ ';
            end
            loglog( ddtt, maxval, ['v-.'], 'MarkerIndices', 1:3:length(ddtt), 'DisplayName', ['$K_w = $',  this, '$\cdot K$ '])
            hold on
            ylabel('$\| \lambda \|$', 'interpreter', 'latex')
            xlabel('t (s)', 'interpreter', 'latex')
            set(gca, 'FontSize', 13)
            drawnow
        end
    end
    
    figure(90+elem)
    
    yy = ylim();
    if (yy(1) > 0.1)
        yy(1) = 0.1;
    end
    ylim(yy);
    
    legend('interpreter', 'latex', 'location', 'best')
    
    drawnow
    yy = ylim();
    yy(1) = 0.1;
    
    xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
    plot(xx, yy, 'k-.', 'HandleVisibility','off')
    ylim(yy);
    drawnow
    print(['ExampleCompressibility-AMatrix-', ElementType], '-dpdf')
    print(['ExampleCompressibility-AMatrix-', ElementType], '-dpng')
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



