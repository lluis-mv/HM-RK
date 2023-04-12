function [] = ExampleThree()
addpath('../Sources')




% 1. Define the problem
indentation = 0.05;

CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 2;



[Nodes, Elements] = ReadTheMesh('Mesh4.msh');

[Nodes1, Elements1] = SimplifyOrder(Nodes,Elements);
nSteps = 100;
dt = 1.0/nSteps;

%MakeSketch
MakeSketch(Nodes1, Elements1);
print('SketchFootingQ', '-dpdf')
clc; clf; close all; 



ind = find(Nodes(:,2) == max( Nodes(:,2)));
xx = sort(Nodes(ind,1));
ind = find(xx == 1);
l = 0.5*(xx(ind)+xx(ind+1));
l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));


if ( true)
tic
[U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodalQuad(Nodes1, Elements1, CP, dt, nSteps, 'Q4Q4', 1);
toc
FF = [information2.F];
FF(1:2:end) = FF(1:2:end)/l;
figure(212); clf
plot( [information2.t]*indentation, FF(1:2:end), 'r', 'linewidth', 2,'DisplayName', ['NS-Q4Q4'])
hold on
figure(214); clf
plot( [information2.t]*indentation, FF(2:2:end), 'r', 'linewidth', 2, 'DisplayName', ['NS-Q4Q4'])
hold on


figure(557); clf
PlotNodal(Nodes1, Elements1, U(3:3:end) )
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
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'Q4Q4', 1);
toc
FF = [information.F];
FF(1:2:end) = FF(1:2:end)/l;
figure(212)
plot( [information.t]*indentation, FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['Q4Q4'])
figure(214)
plot( [information.t]*indentation, FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['Q4Q4'])


figure(556); clf
PlotNodal(Nodes1, Elements1, U(3:3:end) )
drawnow


figure(956); clf
SV = [];
pEff = [];
for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        SV(i,j) = GPInfo(i,j).StressNew(2);
        pEff(i,j) = mean(GPInfo(i,j).StressNew(1:3));
    end
end
PlotHistoryVariable( Nodes1, Elements1, GPInfo, SV);
drawnow



figure(356); clf
PlotHistoryVariable( Nodes1, Elements1, GPInfo, pEff);
drawnow
end

tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'M4Q4', 1);
toc
FF = [information.F];
FF(1:2:end) = FF(1:2:end)/l;
figure(212)
plot( [information.t]*indentation, FF(1:2:end), 'c-.', 'linewidth', 2,'DisplayName', ['Q4Q4Q4'])
figure(214)
plot( [information.t]*indentation, FF(2:2:end), 'c-.', 'linewidth', 2,'DisplayName', ['Q4Q4Q4'])


figure(558); clf
PlotNodal(Nodes1, Elements1, U(4:4:end) )
drawnow


figure(958); clf
SV = [];
pEff = [];
for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        SV(i,j) = GPInfo(i,j).StressNew(2);
        pEff(i,j) = mean(GPInfo(i,j).StressNew(1:3));
    end
end
PlotHistoryVariable( Nodes1, Elements1, GPInfo, SV);
drawnow



figure(358); clf
PlotHistoryVariable( Nodes1, Elements1, GPInfo, pEff);
drawnow





tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'Q8Q4');
toc

FF = [information.F];
FF(1:2:end) = FF(1:2:end)/l2;
figure(212)
plot( [information.t]*indentation, FF(1:2:end), 'b-.', 'linewidth', 2, 'DisplayName',  ['Q8Q4'])
hold on


figure(214)
plot( [information.t]*indentation, FF(2:2:end), 'b-.', 'linewidth', 2, 'DisplayName', ['Q8Q4'])
hold on

figure(559); clf
PlotNodal(Nodes, Elements, U(3:3:end) )
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
for iii = [956:959]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis([-15,-2]);
    colorbar
    drawnow
    pause(1)
    
    fig = figure(iii);
    exportgraphics(fig,['Q1-SV-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(357)
cc = caxis;
i = 1;
pause(1)
for iii = [356:359]
    fig = figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis([-13,-5]);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['Q1-pEff-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(557)
cc = caxis;
i = 1;
pause(1)
for iii = [556:559]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis([0, 22]);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['Q1-Water-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(212)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z$ (m)', 'interpreter', 'latex')
ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
ylim([10,40])
drawnow
fig = figure(212);
exportgraphics(fig,['Q1-Reaction.pdf'], 'BackgroundColor', 'none','ContentType','vector');

figure(214)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
ylim([0,25])
drawnow
xlabel('Footing indentation, $u_z$ (m)', 'interpreter', 'latex')
ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
fig = figure(214);
exportgraphics(fig,['Q1-Water.pdf'], 'BackgroundColor', 'none','ContentType','vector');
