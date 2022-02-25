function [] = CompareThis()

figure(212); hold off; clf;
figure(214); hold off; clf;

figure(1); hold off; clf;
figure(2); hold off; clf;
figure(11); hold off; clf;
figure(12); hold off; clf;


addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = false;
CP.MCC = 2;

eSize = [0.35];

CP.RK = 1;


LinearElastic = false;
if ( LinearElastic)
    CP.Elastic = true;
    CP.MCC = false;
end



model = createpde(1);


R1 = [3,5, 0, 1, 5, 5, 0, 0, 0, 0, -5, -5]';



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


RK = 2;


nSteps = 20;
dt = 0.15/100*nSteps;


%%%%%%%%%%%%%%%%%%%%%%%%% Q8Q4
ElementsQ; NodesQ;


if ( LinearElastic)
    [U1, GPInfo1,  information1] = ComputeImplicitLinearProblem(NodesQQ, ElementsQQ, CP, dt, nSteps, 'Q8Q4',1);
else
    [U1, GPInfo1, rrr1,  information1] = ComputeImplicitNonLinearProblem(NodesQQ, ElementsQQ, CP, dt, nSteps, 'Q8Q4', 0);
end

figure(1)
PlotQ8Nodal(NodesQQ, ElementsQQ, U1(3:3:end) )
drawnow;
colorbar; 

figure(11)
PlotQ8Nodal(NodesQQ, ElementsQQ, U1(2:3:end) )
drawnow;
colorbar;


FF = [information1.F];
figure(212)
plot( [information1.t], FF(1:2:end), 'r', 'linewidth', 2 , 'DisplayName', ['Q8Q4'])
hold on

figure(214)
plot( [information1.t], FF(2:2:end), 'r', 'linewidth', 2, 'DisplayName', ['Q8Q4'])
hold on


%%%%%%%%%%%%%%%%%%%%%%%%% T6T3 




if ( LinearElastic)
    [U2, GPInfo,  information2] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', 1);
else
    [U2, GPInfo, rrr,  information2] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', 0);
end

figure(2)
pdeplot(model,'XYData',U2(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

figure(12)
pdeplot(model,'XYData',U2(2:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;
FF = [information2.F];
figure(212)
plot( [information2.t], FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['T6T3'])
hold on
legend()
figure(214)
plot( [information2.t], FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['T6T3'])
hold on
legend()



return;
%%%%%%%%%%%%%%%%%%%%%%%%% T6T6

if ( LinearElastic)
    [U3, GPInfo2,  information3] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T6', RK, 1);
else
     [U3, GPInfo2, rrr2,  information3] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, 'T6T6', RK, 1, false);
end
figure(3)
pdeplot(model,'XYData',U3(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;



FF = [information3.F];
figure(212)
plot( [information3.t], FF(1:2:end), 'b--', 'linewidth', 2,'DisplayName', ['T6T6'])
hold on

figure(214)
plot( [information3.t], FF(2:2:end), 'b--', 'linewidth', 2, 'DisplayName', ['T6T6'])
hold on

%%%%%%%%%%%%%%%%%%%%%%%%% T6T3  Implicit


[Uimp, GPInfo, rrr,  informationI] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
figure(4)
pdeplot(model,'XYData',Uimp(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

figure(100)
pdeplot(model,'XYData',Uimp(3:3:end)-U2(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

FF = [informationI.F];
figure(212)
plot( [informationI.t], FF(1:2:end), 'k-.', 'linewidth', 2, 'DisplayName',  ['T6T3 Implicit'])
hold on

figure(214)
plot( [informationI.t], FF(2:2:end), 'k-.', 'linewidth', 2, 'DisplayName', ['T6T3 Implicit'])
hold on


% writting everything down
    
maxx = 0;
for i = [1:4]
    figure(i)
    xx = caxis();
    maxx = maxx+xx(2);
end

for i = [1:4]
    figure(i)
    caxis([0, maxx/4])
    axis equal
    xlim([0,4])
    ylim([-4,0])
    print(['WaterPressurePlastic-', num2str(i)], '-dpdf')
end

figure(212)

xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Footing pressure (kN)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
print('Plastic-RR', '-dpdf')

figure(214)
xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Water pressure (kPa)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
print('Plastic-WP', '-dpdf')
