function [] = CompareThis()

figure(212); hold off; clf;
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = true;

eSize = [1.75];


Elem = 2;


if (Elem == 1)
    ElementType = 'T3T3';
elseif (Elem == 2)
    ElementType = 'T6T3';
else
    ElementType = 'T6T6';
end

model = createpde(1);


R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';

g = decsg(R1);
geometryFromEdges(model, g);

if ( Elem == 1)
    mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
else
    mesh = generateMesh(model, 'Hmax', eSize);
end

Nodes = mesh.Nodes';
Elements = mesh.Elements';



% Estimate the element size
modelA = createpde(1);
geometryFromEdges(modelA, g);
mesha = generateMesh(modelA, 'Hmax', eSize, 'GeometricOrder','linear');
Nodesa = mesha.Nodes';
Elementsa = mesha.Elements';
[GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));


RK = 2;

Perme = 10.^(-7);
CP.k = Perme;

nSteps = 30;
dt = 0.15/nSteps;

[U, GPInfo, rrr,  information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1, false);
figure(555)
pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
drawnow;




FF = [information.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), 'DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on



[U2,GPInfo2, rrr2, information2] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType);
figure(556)
pdeplot(model,'XYData',U2(3:3:end),'ColorMap','jet');
drawnow;

figure(557)
pdeplot(model,'XYData', U(3:3:end)- U2(3:3:end),'ColorMap','jet');
drawnow;



FF = [information2.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), 'DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on

