function [] = ExampleThree()
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 2;

eSize= 0.2;

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
[U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information2.F];
figure(212); clf
plot( [information2.t], FF(1:2:end), 'r', 'linewidth', 2,'DisplayName', ['NS-T3T3'])
hold on
figure(214); clf
plot( [information2.t], FF(2:2:end), 'r', 'linewidth', 2, 'DisplayName', ['NS-T3T3'])
hold on


figure(557); clf
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(957); clf
SV = [GPNodes.StressNew];
SV = SV(2,:);
PlotHistoryVariableNodal( Nodes1, Elements1, GPNodes, SV);
drawnow


figure(357); clf
SV = [GPNodes.StressNew];
pEff = mean(SV(1:3,:));
PlotHistoryVariableNodal( Nodes1, Elements1, GPNodes, pEff);
drawnow



tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information.F];
figure(212)
plot( [information.t], FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])
figure(214)
plot( [information.t], FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])


figure(556); clf
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow


figure(956); clf
SV = [GPInfo.StressNew];
SV = SV(2,:)';
PlotHistoryVariable( Nodes1, Elements1, GPInfo, SV);
drawnow



figure(356); clf
SV = [GPInfo.StressNew];
pEff = mean(SV(1:3,:))';
PlotHistoryVariable( Nodes1, Elements1, GPInfo, pEff);
drawnow





tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
toc

FF = [information.F];
figure(212)
plot( [information.t], FF(1:2:end), 'b-.', 'linewidth', 2, 'DisplayName',  ['T6T3'])
hold on


figure(214)
plot( [information.t], FF(2:2:end), 'b-.', 'linewidth', 2, 'DisplayName', ['T6T3'])
hold on

figure(559); clf
pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(959); clf
SV = [];
pEff = [];
for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        SV(i,j) = GPInfo(i,j).StressNew(2);
        pEff(i,j) = mean(GPInfo(i,j).StressNew(1:3));
    end
end
PlotHistoryVariable( Nodes, Elements, GPInfo, SV);
drawnow


figure(359); clf
PlotHistoryVariable( Nodes, Elements, GPInfo, pEff);
drawnow


figure(957)
cc = caxis;
i = 1;
pause(1)
for iii = [956, 957, 959]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)

    print(['F1-SV-', num2str(i)], '-dpdf');
    i = i+1;
end


figure(357)
cc = caxis; 
i = 1;
pause(1)
for iii = [356, 357, 359]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    print(['F1-pEff-', num2str(i)], '-dpdf');
    i = i+1;
end


figure(557)
cc = caxis; 
i = 1;
pause(1)
for iii = [556, 557, 559]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    print(['F1-Water-', num2str(i)], '-dpdf');
    i = i+1;
end


figure(212)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
print('F1-Reaction', '-dpdf')

figure(214)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
print('F1-Water', '-dpdf')