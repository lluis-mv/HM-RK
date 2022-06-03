function [] = ExampleCompressibility()
close all
addpath('../Sources')
% 1. Define the problem




CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 0.01;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

CP.Compressibility = 0;



for i = 1:2
    figure(i);
    clf;
    newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    colororder(newcolors)
end

[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);



ElementType = 'T6T3';
Nodes = NodesT;
Elements = ElementsT;
ThisNumber = 6;

nSteps = 100;
t = 1;


K = CP.E/3/(1-2*nu);
for C = [10, 100, 1000, 10000]*K
    CP.Compressibility = 1/C;
    [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, t/nSteps, nSteps, ElementType, 1, 0);
    
    
    figure(1)
    index = find(Nodes(:,1) == 0);
    Uv = U(3*(index-1)+2);
    y = Nodes(index,2);
    [ya, index] = sort(y);
    Uv = Uv(index);
    % [Uv, ya] = CorrectInterpolation(Uv, ya);
    plot(Uv, ya, '', 'linewidth', 2, 'DisplayName', 'Reference')
    hold on
    drawnow
    
    
    figure(2)
    index = find(Nodes(:,1) == 0);
    WP = U(3*(index-1)+3);
    y = Nodes(index,2);
    [ya, index] = sort(y);
    WPa = WP(index);
    % [WPa, ya] = CorrectInterpolation(WPa, ya);
    plot(WPa, ya, '', 'linewidth', 2, 'DisplayName', 'Reference')
    hold on
    drawnow
end
