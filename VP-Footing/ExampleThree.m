function [] = ExampleThree()
addpath('../Sources')


clear all; clc; clf; close all;


% 1. Define the problem
indentation = 0.1/3600;

CP.HydroMechanical = true;
CP.E = 1000;
% CP.nu = 0.3;
% nu = CP.nu;
% CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 3;

CP.kappa = 0.01;
CP.lambda = 0.1;
CP.M_MCC = 1;
CP.nu = 0.3;

CP.n = 3;
CP.r = 4;

CP.m = 1.75;

CP.PerzynaN = 1;
CP.PerzynaEta = 1000;
CP.RK = -6;

eSize= 0.15;
MakeSketch(eSize)
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



nSteps = 400;
dt = 3600.0/nSteps;




ind = find(Nodes(:,2) == max( Nodes(:,2)));
xx = sort(Nodes(ind,1));
ind = find(xx == 1);
l = 0.5*(xx(ind)+xx(ind+1));
l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));


iCase = 1;

COLOR = ['krgbcm']

for ETA = [0, 100, 500, 1000]


    color = COLOR(iCase);
    CP.PerzynaEta = ETA;
    CP.MCC = 3;
    if ( CP.PerzynaEta == 0)
        CP.MCC = 4;
    end



    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
    toc

    FF = [information.F];
    FF(1:2:end) = FF(1:2:end)/l2;
    figure(212)
    plot( [information.t]*indentation, FF(1:2:end), color, 'linewidth', 2, 'DisplayName',  ['$\eta = $', num2str(ETA)])
    hold on


    figure(214)
    plot( [information.t]*indentation, FF(2:2:end), color, 'linewidth', 2, 'DisplayName',  ['$\eta = $', num2str(ETA)])
    hold on

    figure(500+iCase); clf
    pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
    drawnow

    figure(900+iCase); clf
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


    figure(300+iCase); clf
    PlotHistoryVariable( Nodes, Elements, GPInfo, pEff);
    drawnow
    GenerateFigures(iCase)

    iCase = iCase+1;

end

function [] = GenerateFigures(iCase)

figure(901)
cc = caxis;
cc(2) = 0.0;
i = 1;
pause(1)
for iii = [901:(900+iCase)]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)

    fig = figure(iii);
    exportgraphics(fig,['F1-SV-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(301)
cc = caxis;
cc(2) = 0.0;
i = 1;
pause(1)
for iii = [301:(300+iCase)]
    fig = figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['F1-pEff-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(501)
cc = caxis;
i = 1;
pause(1)
for iii = [501:(500+iCase)]
    figure(iii)
    axis equal; xlim([0,4]); ylim([-4, 0]); axis off
    colormap jet
    caxis(cc);
    colorbar
    drawnow
    pause(1)
    fig = figure(iii);
    exportgraphics(fig,['F1-Water-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
    i = i+1;
end


figure(212)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z$ (m)', 'interpreter', 'latex')
ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
fig = figure(212);
exportgraphics(fig,['F1-Reaction.pdf'], 'BackgroundColor', 'none','ContentType','vector');

figure(214)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
xlabel('Footing indentation, $u_z$ (m)', 'interpreter', 'latex')
ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
fig = figure(214);
exportgraphics(fig,['F1-Water.pdf'], 'BackgroundColor', 'none','ContentType','vector');
