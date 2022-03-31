function TestThis()
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = false;
CP.MCC = false;

ESIZE = 0.35;
eSize = ESIZE;
RKReference = 8;
RKMethods = [8,1:7];

RKReference = 8;


Elem = 1;


if (Elem == 1)
    ElementType = 'T3T3';
     ThisNumber = 200;
elseif (Elem == 2)
    ElementType = 'T6T3';
else
    ElementType = 'T6T6';
end

model = createpde(1);


R1 = [3,5, 0, 1, 3, 3, 0, 0, 0, 0, -3, -3]';

g = decsg(R1);
geometryFromEdges(model, g);

if ( Elem == 1)
    mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
else
    mesh = generateMesh(model, 'Hmax', eSize);
end

Nodes = mesh.Nodes';
Elements = mesh.Elements';

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));


ddtt = 10.^linspace(0,2,40);

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
    
    
    
end
figure(j+1);
clf;

loglog(ddtt, minval2, 'm*-.', ddtt, maxval2, 'c*-.')

drawnow;
hold on
loglog(ddtt, minval, 'r*-.', ddtt, maxval, 'b*-.')
drawnow

xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
ylabel('$\| \lambda \|$', 'interpreter', 'latex')
set(gca, 'FontSize', 14)
drawnow
yy = ylim();
xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
plot(xx, yy, 'k-.')
ylim(yy);
ll = legend('min$(|\lambda|)$ Primal', 'max$(|\lambda|)$ Primal', ...
    'min$(|\lambda|)$ Stab', 'max$(|\lambda|)$ Stab', 'location', 'best');
set(ll, 'interpreter', 'latex')
xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])

print(['ExampleOne-Radii-', ElementType], '-dpdf')




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