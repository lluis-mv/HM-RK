function [] = ExampleFour()
addpath('../Sources')


clear all; clc; clf; close all;


% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.k = 1E-12;
CP.k = 1E-7;
CP.Elastic = false;
CP.MCC = 3;


CP.kappa = 0.01;
CP.lambda = 0.1;
CP.M_MCC = 1.3;
CP.nu = 0.3;

CP.n = 1.5;
CP.r = 2;

CP.m = 1.75;

CP.PerzynaN = 1;
CP.PerzynaEta = 1000;
CP.RK = -2;

[Nodes, Elements] = ReadTheMesh('Mesh.msh');



dt = 1/25;
nSteps = 300;



% ind = find(Nodes(:,2) == max( Nodes(:,2)));
% xx = sort(Nodes(ind,1));
% ind = find(xx == 1);
% l = 0.5*(xx(ind)+xx(ind+1));
% l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));
l2 = 1;

iCase = 1;

COLOR = ['krgbcm']

for ETA = [0, 10, 100, 1000]


    color = COLOR(iCase);
    col = COLOR(iCase);
    if ( color == 'r')
        color = 'r-.';
    end
    CP.PerzynaEta = ETA;
    CP.MCC = 3;
    if ( CP.PerzynaEta == 0)
        CP.MCC = 4;
    end
    
    


    tic
    [U, GPInfo, rrr,  information, ~, ErrorNorms] = ComputeImplicitNonLinearProblemContact(Nodes, Elements, CP, dt, nSteps, 'T6T3');
    toc

   NameToDisplay = ['$\eta = $', num2str(ETA), ' s'];
    if ( ETA == 0)
        NameToDisplay = 'Elasto-Plastic';
    end

    FF = [information.F];
    FF(1:2:end) = FF(1:2:end)/l2;
    figure(212)
    semilogx( [information.t], FF(1:2:end), color, 'linewidth', 2, 'DisplayName',  NameToDisplay)
    hold on


    figure(213)
    semilogx( [information.t], [information.u], color, 'linewidth', 2, 'DisplayName',  NameToDisplay)
    hold on


    figure(214)
    semilogx( [information.t], FF(2:2:end), color, 'linewidth', 2, 'DisplayName',  NameToDisplay)
    hold on

%     figure(500+iCase); clf
%     pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
%     drawnow
% 
%     figure(900+iCase); clf
%     SV = [];
%     pEff = [];
%     for i = 1:size(GPInfo,1)
%         for j = 1:size(GPInfo, 2)
%             SV(i,j) = -GPInfo(i,j).StressNew(2);
%             pEff(i,j) = -mean(GPInfo(i,j).StressNew(1:3));
%         end
%     end
%     PlotHistoryVariable( Nodes, Elements, GPInfo, SV);
%     drawnow
% 
% 
%     figure(300+iCase); clf
%     PlotHistoryVariable( Nodes, Elements, GPInfo, pEff);
%     drawnow
    GenerateFigures(iCase)


%     figure(20);
%     nnn = 0:(length(ErrorNorms.Iter0)-1);
%     semilogy( nnn, ErrorNorms.Iter0, [col, '*-.'], 'linewidth', 2, ...
%         'DisplayName',  ['$\eta = $', num2str(ETA)]);
%     hold on
%     xlabel('Iteration, $i$', 'interpreter', 'latex')
%     ylabel('Norm of the residual, $\|\mathbf{R}_i\|$', 'interpreter', 'latex')
%     set(gca, 'FontSize', 15)
%     legend('interpreter', 'latex', 'location', 'best')
%     drawnow
%     print(['Contact-C-Res0-'], '-dpdf')
% 
% 
%     fig = figure(21);
%     loglog( ErrorNorms.Iter0(1:end-1), ErrorNorms.Iter0(2:end), [col, '*-.'], 'linewidth', 2, ...
%         'DisplayName',  ['$\eta = $', num2str(ETA)]);
%     hold on
%     xlabel('$$\|\mathbf{R}_i\|$', 'interpreter', 'latex')
%     ylabel('$$\|\mathbf{R}_{i+1}\|$', 'interpreter', 'latex')
%     set(gca, 'FontSize', 15)
%     legend('interpreter', 'latex', 'location', 'best','location', 'northwest')
%     drawnow
%     exportgraphics(fig,'Contact-C-Res01-.pdf', 'BackgroundColor', 'none','ContentType','vector');
% 
% 
    figure(30);
    nnn = 0:(length(ErrorNorms.IterEnd)-1);
    semilogy( nnn, ErrorNorms.IterEnd, [col, '*-.'], 'linewidth', 2, ...
        'DisplayName',  ['$\eta = $', num2str(ETA)]);
    hold on
    xlabel('Iteration', 'interpreter', 'latex')
    ylabel('NormResidual', 'interpreter', 'latex')
    set(gca, 'FontSize', 15)
    legend('interpreter', 'latex', 'location', 'best')
    drawnow
    print(['Contact-C-ResF-'], '-dpdf')

    fig = figure(31);
    loglog( ErrorNorms.IterEnd(1:end-1), ErrorNorms.IterEnd(2:end), [col, '*-.'], 'linewidth', 2, ...
        'DisplayName',  ['$\eta = $', num2str(ETA)]);
    hold on
    xlabel('$R_i$', 'interpreter', 'latex')
    ylabel('$R_{i+1}$', 'interpreter', 'latex')
    set(gca, 'FontSize', 15)
    legend('interpreter', 'latex', 'location', 'best','location', 'northwest')
    drawnow
    exportgraphics(fig,'Contact-C-ResF1-.pdf', 'BackgroundColor', 'none','ContentType','vector');




    iCase = iCase+1;

end

% fig = figure(21);
% xx = xlim();
% yy = ylim();
% DrawLines();
% xlim(xx);
% ylim(yy);
% drawnow
% exportgraphics(fig,'Contact-C-Res01-.pdf', 'BackgroundColor', 'none','ContentType','vector');
% 
fig = figure(31);
xx = xlim();
yy = ylim();
DrawLines();
xlim(xx);
ylim(yy);
drawnow
exportgraphics(fig,'Contact-C-ResF1-.pdf', 'BackgroundColor', 'none','ContentType','vector');

function [] = GenerateFigures(iCase)

% figure(901)
% cc = caxis;
% i = 1;
% pause(1)
% for iii = [901:(900+iCase)]
%     figure(iii)
%     axis equal; %xlim([0,4]); ylim([-4, 0]); 
%     axis off
%     colormap jet
%     caxis(cc);
%     colorbar('southoutside')
%     drawnow
%     pause(1)
% 
%     fig = figure(iii);
%     exportgraphics(fig,['C1-SV-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
%     i = i+1;
% end
% 
% 
% figure(301)
% cc = caxis;
% i = 1;
% pause(1)
% for iii = [301:(300+iCase)]
%     fig = figure(iii)
%     axis equal; %xlim([0,4]); ylim([-4, 0]); 
%     axis off
%     colormap jet
%     caxis(cc);
%     colorbar('southoutside')
%     drawnow
%     pause(1)
%     fig = figure(iii);
%     exportgraphics(fig,['C1-pEff-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
%     i = i+1;
% end
% 
% 
% figure(501)
% cc = caxis;
% i = 1;
% pause(1)
% for iii = [501:(500+iCase)]
%     figure(iii)
%     axis equal; %xlim([0,4]); ylim([-4, 0]); 
%     axis off
%     colormap jet
%     caxis(cc);
%     colorbar('southoutside')
%     drawnow
%     pause(1)
%     fig = figure(iii);
%     exportgraphics(fig,['C1-Water-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
%     i = i+1;
% end


figure(212)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
legend('location', 'south')
xlabel('Time, $t$ (s)', 'interpreter', 'latex')
ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
fig = figure(212);
exportgraphics(fig,['C1-Reaction.pdf'], 'BackgroundColor', 'none','ContentType','vector');

figure(213)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
legend('location', 'south')
xlabel('Time, $t$ (s)', 'interpreter', 'latex')
ylabel('Vertical displacement, $u_z$ (m)', 'interpreter', 'latex')
fig = figure(213);
exportgraphics(fig,['C1-Displacement.pdf'], 'BackgroundColor', 'none','ContentType','vector');


figure(214)
legend('location', 'best', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
legend('location', 'south')
xlabel('Time, $t$ (s)', 'interpreter', 'latex')
ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
fig = figure(214);
exportgraphics(fig,['C1-Water.pdf'], 'BackgroundColor', 'none','ContentType','vector');


function [] = DrawLines()

xx = 10.^linspace(-20,20,5);
for c = (-20:2:20)
    yy = 10.^(2*log10(xx)+c);
    zz = 10.^(-10) * (1:length(yy));
    plot3(xx,yy,  zz, 'k-.', 'HandleVisibility', 'off')
    hold on
end
