function [] = ExampleOne()

close all
addpath('../Sources')
% 1. Define the problem

T = 1E-7;


CP.HydroMechanical = true;
CP.E = 100;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.04;

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
print('ExampleOne-FemMesh', '-dpdf')

Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');


NStepsRef = 1;



Stab = 1;
ElementType = 'T6T3';
Nodes = Nodes2;
Elements = Elements2;
Color = 'b';
dt = 1;

[U, GPInfo] = ComputeImplicitNonLinearProblemRoot(Nodes, Elements, CP, dt, NStepsRef, ElementType, Stab);
figure(2)
PlotNodal(Nodes, Elements, U(2:3:end))
axis equal
axis off
colorbar



figure(3); clf
SV = [];

for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        SV(i,j) = GPInfo(i,j).StrainNew(2);

    end
end
PlotHistoryVariable( Nodes, Elements, GPInfo, SV);
drawnow


axis equal
axis off
colorbar
colormap jet


for i = 1:3; figure(i); drawnow; pause(1); print(['Figure', num2str(i)], '-dpng'); end