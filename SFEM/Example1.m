function [] = ExampleOne()

addpath('../')
% 1. Define the problem

T = 1E-02;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;

eSize = 0.045;
% eSize = 10.0;

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



NStepsRef = 10;
dt = t/NStepsRef;

[U1,GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt,NStepsRef, 'T3T3');
index = find(Nodes(:,1) == 0);
WP = U1(3*(index-1)+3);
y = Nodes(index,2);
[y, index] = sort(y);
WP = WP(index);
figure(3); clf; 
plot(WP, y, 'k', 'linewidth', 2, 'DisplayName', 'T3T3')
hold on

[U2,GPInfo] = ComputeImplicitNonLinearProblem(Nodes2, Elements2, CP, dt,NStepsRef, 'T6T3');
index = find(Nodes(:,1) == 0);
WP = U2(3*(index-1)+3);
y = Nodes(index,2);
[y, index] = sort(y);
WP = WP(index);
plot(WP, y, 'r-.', 'linewidth', 2, 'DisplayName', 'T6T3')

[U3,GPInfo] = ComputeImplicitNonLinearProblemNodal(Nodes, Elements, CP, dt,NStepsRef, 'T3T3');
index = find(Nodes(:,1) == 0);
WP = U3(3*(index-1)+3);
y = Nodes(index,2);
[y, index] = sort(y);
WP = WP(index);
plot(WP, y, 'g-.', 'linewidth', 2, 'DisplayName', 'T6T3')


