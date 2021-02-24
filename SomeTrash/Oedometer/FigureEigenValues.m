function [] = FigureEigenValues()

addpath('../')
% 1. Define the problem



CP.E = 1000;
CP.nu = 0.3;
CP.k = 1E-4;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




eSize = 0.05;

model = createpde(1);

dx = 0.1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes1 = mesh.Nodes';
Elements1 = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';




NSteps = 10.^linspace(0, 5, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);

i = 1;

firstTime = true;

ddtt = 10.^linspace(-7,1, 10);


for j = 1:3
    figure(j);
    clf;
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
        loglog(ddtt, minval, 'r*-.', ddtt, maxval, 'b*-.')
        drawnow;
        
        
        
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
        
        hold on
        semilogx(ddtt, minval2, 'm*-.', ddtt, maxval2, 'c*-.')
        drawnow
        hold off
    end
    
end


function [A, nDirichlet, nNoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, RKMethod, AlphaStabM)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);


[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);


[C, K, X, f, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, C, K);


PostProcessResults(Nodes, Elements, X, GPInfo, 0, true, ['ThisProblem-', ElementType]);


A = C\(K);

nNoDir = [];
for i = 1:3*nNodes
    if ( any(i == nDirichlet))
        % do nothing
    else
        nNoDir = [nNoDir, i];
    end
end







