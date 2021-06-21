function [] = ExampleThree()
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-2;
CP.Elastic = false;
CP.MCC = 2;

eSize= 0.35;



model = createpde(1);


R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';



g = decsg(R1);
geometryFromEdges(model, g);
mesh = generateMesh(model, 'Hmax', eSize);
Nodes = mesh.Nodes';
Elements = mesh.Elements';


model1 = createpde(1);
geometryFromEdges(model1, g);
mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
Nodes1 = mesh1.Nodes';
Elements1 = mesh1.Elements';



nSteps = 100;
dt = 0.15/nSteps;

tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information.F];
figure(212); clf;
plot( [information.t], FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])
hold on

figure(214); clf;
plot( [information.t], FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])
hold on




figure(556)
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow


figure(956)
SV = [GPInfo.StressNew];
SV = SV(2,:);
PlotHistoryVariable( Nodes1, Elements1, GPInfo, SV);
drawnow


tic

[U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information2.F];
figure(212)
plot( [information2.t], FF(1:2:end), 'b--', 'linewidth', 2,'DisplayName', ['S-T3T3'])
hold on
legend('location', 'best')
figure(214)
plot( [information2.t], FF(2:2:end), 'b--', 'linewidth', 2, 'DisplayName', ['S-T3T3'])
legend('location', 'best')
hold on


figure(557)
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(957)
SV = [GPNodes.StressNew];
SV = SV(2,:);
PlotHistoryVariableNodal( Nodes1, Elements1, GPNodes, SV);
drawnow


figure(958)
pdeplot(model1,'XYData',SV,'ColorMap','jet');
drawnow





tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
toc

FF = [information.F];
figure(212)
plot( [information.t], FF(1:2:end), 'k-.', 'linewidth', 2, 'DisplayName',  ['T6T3'])
hold on


figure(214)
plot( [information.t], FF(2:2:end), 'k-.', 'linewidth', 2, 'DisplayName', ['T6T3'])
hold on

figure(559)
pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
drawnow



figure(956)
m2 = caxis();

figure(957)
caxis(m2)


