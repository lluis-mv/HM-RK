function [] = CompareElementsElastic()

figure(212); hold off; clf;
figure(214); hold off; clf;
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = true;
CP.MCC = true;

eSize = [0.35];




LinearElastic = false;
if ( LinearElastic)
    CP.Elastic = true;
    CP.MCC = false;
end



[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);

RK = 2;


nSteps = 100;
dt = 0.15/nSteps;


%%%%%%%%%%%%%%%%%%%%%%%%% T3T3 

[U1, GPInfo1, rrr1,  information1] = ComputeNLProblem(NodesT, ElementsT, CP, dt, nSteps, 'T6T3', RK, 1, false);



figure(1)
PlotNodal(NodesT, ElementsT, U1(3:3:end))
drawnow;
colorbar; 




FF = [information1.F];
figure(212)
plot( [information1.t], FF(1:2:end), 'r', 'linewidth', 2 , 'DisplayName', ['T6T3'])
hold on
drawnow

figure(214)
plot( [information1.t], FF(2:2:end), 'r', 'linewidth', 2, 'DisplayName', ['T6T3'])
hold on
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%% Q8Q4





[U2, GPInfo, rrr,  information2] = ComputeNLProblem(NodesQ, ElementsQ, CP, dt, nSteps, 'Q8Q4', RK, 1, false);


figure(2)
PlotNodal(NodesQ, ElementsQ, U2(3:3:end))
drawnow;
colorbar; 

FF = [information2.F];
figure(212)
plot( [information2.t], FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['Q8Q4'])
hold on
drawnow

figure(214)
plot( [information2.t], FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['Q8Q4'])
hold on
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%% T6T3  Implicit


[Uimp, GPInfo, rrr,  informationI] = ComputeImplicitNonLinearProblem(NodesT, ElementsT, CP, dt, nSteps, 'T6T3');

figure(3)
PlotNodal(NodesT, ElementsT, Uimp(3:3:end))
drawnow;
colorbar; 


FF = [informationI.F];
figure(212)
plot( [informationI.t], FF(1:2:end), 'k:', 'linewidth', 2, 'DisplayName',  ['T6T3 Implicit'])
hold on
drawnow

figure(214)
plot( [informationI.t], FF(2:2:end), 'k:', 'linewidth', 2, 'DisplayName', ['T6T3 Implicit'])
hold on
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%% T6T3  Implicit


[Uimp2, GPInfo, rrr,  informationI] = ComputeImplicitNonLinearProblem(NodesQ, ElementsQ, CP, dt, nSteps, 'Q8Q4');

figure(4)
PlotNodal(NodesQ, ElementsQ, Uimp2(3:3:end) )
drawnow;
colorbar; 


FF = [informationI.F];
figure(212)
plot( [informationI.t], FF(1:2:end), 'b:', 'linewidth', 2, 'DisplayName',  ['Q8Q4 Implicit'])
hold on
drawnow

figure(214)
plot( [informationI.t], FF(2:2:end), 'b:', 'linewidth', 2, 'DisplayName', ['Q8Q4 Implicit'])
hold on
drawnow






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
    axis off
    pbaspect([1 1 10])
    print(['WaterPressureElastic-', num2str(i)], '-dpdf')
end

figure(212)

xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Footing pressure (kN)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
print('Elastic-RR', '-dpdf')

figure(214)
xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Water pressure (kPa)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
print('Elastic-WP', '-dpdf')
