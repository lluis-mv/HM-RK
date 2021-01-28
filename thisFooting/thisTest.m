function [] = ExampleThree()

figure(30); clf;
figure(50); clf;
figure(900); clf;

addpath('../')
% 1. Define the problem

T = 1E-1;

CP.HydroMechanical = true;
CP.E = 100;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-3;




eSize = [0.15];


figure(50); clf;
RKMethod = 1;
Elem = 2;

% for Elem = [1, 2, 3]
RKReference = 8;
RKMethods = [RKReference, 1,2,4, 6, -1];
RKReference = 8;
RKMethods = [8,1,2,4,6];


Elem = 1;



if (Elem == 1)
    ElementType = 'T3T3';
    ThisNumber = 200;
elseif (Elem == 2)
    ElementType = 'T6T3';
    ThisNumber = 6;
else
    ElementType = 'T6T6';
    ThisNumber = 2000;
end

dx = 0.3; dy = 1;
model = createpde(1);


R1 = [3,5, 0, 1, 3, 3, 0, 0, 0, 0, -3, -3]';
%         R1 = [3,4, 0, 1, 1, 0,  0, 0, -3, -3]';
g = decsg(R1);
geometryFromEdges(model, g);

if ( Elem == 1)
    mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
else
    mesh = generateMesh(model, 'Hmax', eSize);
end

Nodes = mesh.Nodes';
Elements = mesh.Elements';


% First part. compute the eigenvalues
figure(1);
clf;
if ( Elem == 1)
    triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
end
drawnow
axis equal
axis off


% Estimate the element size

mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
Nodesa = mesha.Nodes';
Elementsa = mesha.Elements';
[GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));










nSteps = 10;
RK = 3;

dt = 0.1/nSteps;
[U2,GPInfo2,r2,  information2] = ComputeNLProblem2(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1);
[U,GPInfo,r1,  information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1);









