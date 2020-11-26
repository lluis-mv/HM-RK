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



