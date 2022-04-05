function [] = ExampleOne()
close all
addpath('../')
% 1. Define the problem

T = 0.05;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;



[NodesQ, ElementsQ] = ReadTheMesh('MeshUglyQ.msh');
[NodesT, ElementsT] = ReadTheMesh('MeshUglyT.msh');

[GPInfo] = ComputeElementalMatrices(NodesT, ElementsT, CP, 'T6T3');



NSteps = 10.^linspace(0, 6, 20);
NSteps = floor(NSteps); NSteps = sort(NSteps);
NStepsRef = 1;
ddtt = t./NSteps;

ticks = [min(ddtt), max(ddtt)];
ticks = floor(log10(ticks));
ticks = 10.^(ticks(1):2:ticks(end));
ticks = [ticks, min(ddtt), max(ddtt)];
ticks = 10.^unique( log10(ticks));
ticks = unique(ticks);


for Mesh = [2,3]

    if ( Mesh == 1)
        [NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
        [NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);
    elseif (Mesh == 2)
        [NodesQ, ElementsQ] = ReadTheMesh('ThisMeshQ.msh');
        [NodesT, ElementsT] = ReadTheMesh('ThisMeshT.msh');
    elseif (Mesh == 3)
        [NodesQ, ElementsQ] = ReadTheMesh('MeshUglyQ.msh');
        [NodesT, ElementsT] = ReadTheMesh('MeshUglyT.msh');
    end
    
    
    % Now lets check the RK methods
    for j = [2, 1]
        
        if ( j == 1)
            ElementType = 'T6T3';
            Nodes = NodesT;
            Elements = ElementsT;
        elseif (j == 2)
            ElementType = 'Q8Q4';
            Nodes = NodesQ;
            Elements = ElementsQ;
        end
        
        Stab = 1;
        
        
        for RK = [1,3,8]
            firstTime = true;
            i = 1;
            
            for nSteps = NSteps
                
                dt = t/nSteps;
                [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, Stab);
                if ( firstTime)
                    [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
                end
                [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo);
                
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

            loglog( ddtt, L2, [color, '*-.'], 'DisplayName',  ['$\| p_w-p_w^h\|_{L2}$ RK-', num2str(RK)]);
            hold on
            loglog( ddtt, L2U, [color, 'v-.'], 'DisplayName',  ['$\| \mathbf{u}-\mathbf{u}_h\|_{L2}$ RK-', num2str(RK)]);
            hold on
            xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            print(['ExampleOneUgly-RK-', ElementType, '-MESH-', num2str(Mesh)], '-dpdf')
            
        end
        figure(20+3*j)
        drawnow
        xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
        xticks(ticks);
        legend('location', 'best', 'interpreter', 'latex')
        print(['ExampleOneUgly-RK-', ElementType, '-MESH-', num2str(Mesh)], '-dpdf')
        
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
