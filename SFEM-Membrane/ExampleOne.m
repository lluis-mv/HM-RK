function [] = ExampleOne()

close all
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.0;
CP.k = 1E-8;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




eSize = 10;

model = createpde(1);



R1 = [3,4,0, 48, 48, 0, 0, 44, 44+16, 44]';
% R1 = [3,4,0, 48, 48, 0, 0, 0, 20, 20]';



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


% First part. compute the eigenvalues
figure(1);
clf;
triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k');
drawnow
axis equal
axis off
print('ExampleMembrane-FemMesh', '-dpdf')



% Estimate the element size






dt = 0.01;
nSteps = 1;


tic
[U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information2.F];
FF(1:2:end) = FF(1:2:end);
figure(212); clf
plot( [information2.t], FF(1:2:end), 'r', 'linewidth', 2,'DisplayName', ['NS-T3T3'])
hold on
figure(214); clf
plot( [information2.t], FF(2:2:end), 'r', 'linewidth', 2, 'DisplayName', ['NS-T3T3'])
hold on


figure(557); clf
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(857); clf
pdeplot(model1,'XYData',U(2:3:end),'ColorMap','jet');
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

figure(157); clf;
triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k');
hold on
quiver(Nodes1(:,1), Nodes1(:,2), U(1:3:end), U(2:3:end), 'k')


tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
toc
FF = [information.F];
FF(1:2:end) = FF(1:2:end);
figure(212)
plot( [information.t], FF(1:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])
figure(214)
plot( [information.t], FF(2:2:end), 'g', 'linewidth', 2,'DisplayName', ['T3T3'])


figure(556); clf
pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(856); clf
pdeplot(model1,'XYData',U(2:3:end),'ColorMap','jet');
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

figure(156); clf;
triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k');
hold on
quiver(Nodes1(:,1), Nodes1(:,2), U(1:3:end), U(2:3:end), 'k')



tic
[U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
toc

FF = [information.F];
FF(1:2:end) = FF(1:2:end);
figure(212)
plot( [information.t], FF(1:2:end), 'b-.', 'linewidth', 2, 'DisplayName',  ['T6T3'])
hold on


figure(214)
plot( [information.t], FF(2:2:end), 'b-.', 'linewidth', 2, 'DisplayName', ['T6T3'])
hold on

figure(559); clf
pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
drawnow

figure(859); clf
pdeplot(model,'XYData',U(2:3:end),'ColorMap','jet');
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

figure(159); clf;
triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k');
hold on
quiver(Nodes(:,1), Nodes(:,2), U(1:3:end), U(2:3:end), 'k')



figure(359); clf
PlotHistoryVariable( Nodes, Elements, GPInfo, pEff);
drawnow


figure(957)
cc = caxis;
i = 1;
pause(1)
for iii = [956, 957, 959]
    figure(iii)
    axis equal; xlim([0,48]); ylim([0, 44+16]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    
    fig = figure(iii);
    exportgraphics(fig,['F1-SV-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(357)
cc = caxis;
i = 1;
pause(1)
for iii = [356, 357, 359]
    fig = figure(iii)
    axis equal; xlim([0,48]); ylim([0, 44+16]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['F1-pEff-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(557)
cc = caxis;
i = 1;
pause(1)
for iii = [556, 557, 559]
    figure(iii)
    axis equal; xlim([0,48]); ylim([0, 44+16]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['F1-Water-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(857)
cc = caxis;
i = 1;
pause(1)
for iii = [856, 857, 859]
    figure(iii)
    axis equal; xlim([0,48]); ylim([0, 44+16]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['F1-UY-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end

figure(212)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
fig = figure(212);
exportgraphics(fig,['F1-Reaction.pdf'], 'BackgroundColor', 'none','ContentType','vector');

figure(214)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
fig = figure(214);
exportgraphics(fig,['F1-Water.pdf'], 'BackgroundColor', 'none','ContentType','vector');



function [p, y] = CorrectInterpolation(p1, y1)

alfa = linspace(0, 1);
y = [];
p = [];

for ind = 1:2:length(y1)-1
    pp = p1(ind:ind+2);
    yy = y1(ind:ind+2);
    
    N = [(1 - alfa).*(1-2*alfa);
        4*(1-alfa).*alfa;
        alfa.*(2*alfa-1)];
    y = [y; N'*yy];
    p = [p; N'*pp];
end



function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% Other analytical solution...
nNodes = size(Nodes,1);
M = CP.M;
k = CP.k;
for nod = 1:nNodes
    xx = 1-Nodes(nod,2);
    TT = M * t*k;
    pw = 0;
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw = pw + 2/aux * sin( aux * xx) * exp( - aux^2 * TT);
    end
    Xa(3*(nod-1)+3) = pw;
end


