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


MakeSketch('Mesh-Coarse.msh')
% model = createpde(1);


% R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';



% g = decsg(R1);
% geometryFromEdges(model, g);
% mesh = generateMesh(model, 'Hmax', eSize);
% Nodes = mesh.Nodes';
% Elements = mesh.Elements';

[Nodes, Elements] = ReadTheMesh('Mesh-Coarse.msh');

nSteps = 200;
dt = 3600.0/nSteps;
nSteps = 50;




ind = find(Nodes(:,2) == max( Nodes(:,2)));
xx = sort(Nodes(ind,1));
ind = find(xx == 1);

l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));




COLOR = ['krgbcm']

xAxis = [1,2,4,6,8];

ETAS = [0, 1, 500, 1000];
Resistance = nan*ones(length(xAxis), length(ETAS));
TIME = nan*ones(length(xAxis), length(ETAS));
WaterP = nan*ones(length(xAxis), length(ETAS));

jCase = 1;
for RK = -xAxis
    iCase = 1;

    figure(20); clf;
    figure(21); clf;

    figure(30); clf;
    figure(31); clf;

    for ETA = ETAS


        CP.RK = RK;
        CP.PerzynaEta = ETA;
        CP.MCC = 3;
        if ( CP.PerzynaEta == 0)
            CP.MCC = 4;
        end


        col = COLOR(iCase);
        tic
        [U, GPInfo, rrr,  information, ~, ErrorNorms] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
        TIME(jCase, iCase) = toc;

        FF = [information.F];
        
        Resistance(jCase, iCase) = FF(end-1)/l2;
        WaterP(jCase, iCase) = FF(end);
        

        figure(40); clf

        for iii = 1:4
        plot(xAxis, Resistance(:,iii), [COLOR(iii), '*-.'], 'linewidth', 2);
        hold on
        end
        xlabel('RK order', 'interpreter', 'latex')
        ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        legend('Elasto-Plastic', '$\eta = 1$ s',  '$\eta = 500$ s',  '$\eta = 1000$ s', 'interpreter', 'latex',...
            'location', 'best')
        drawnow
        print('OrderCost-1', '-dpdf')

        figure(41); clf
        for iii = 1:4
        plot(xAxis, WaterP(:,iii), [COLOR(iii), '*-.'], 'linewidth', 2);
        hold on
        end
        xlabel('RK order', 'interpreter', 'latex')
        ylabel('Water pressure (kPa)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        legend('Elasto-Plastic', '$\eta = 1$ s',  '$\eta = 500$ s',  '$\eta = 1000$ s', 'interpreter', 'latex',...
            'location', 'best')
        drawnow
        print('OrderCost-2', '-dpdf')

        figure(42); clf;
        for iii = 1:4
        plot(xAxis, TIME(:,iii), [COLOR(iii), '*-.'], 'linewidth', 2);
        hold on
        end
        xlabel('RK order', 'interpreter', 'latex')
        ylabel('Computational cost (s)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        legend('Elasto-Plastic', '$\eta = 1$ s',  '$\eta = 500$ s',  '$\eta = 1000$ s', 'interpreter', 'latex',...
            'location', 'best')
        drawnow
        print('OrderCost-3', '-dpdf')

%         figure(20); 
%         semilogy( ErrorNorms.Iter0, [col, '*-.'], 'linewidth', 2);
%         hold on
%         xlabel('Iteration', 'interpreter', 'latex')
%         ylabel('NormResidual', 'interpreter', 'latex')
%         set(gca, 'FontSize', 15)
%         print(['OrderCost-Res0-', num2str(jCase)], '-dpdf')
% 
%         figure(21);
%         loglog( ErrorNorms.Iter0(1:end-1), ErrorNorms.Iter0(2:end), [col, '*-.'], 'linewidth', 2);
%         hold on
%         xlabel('$R_i$', 'interpreter', 'latex')
%         ylabel('$R_{i+1}$', 'interpreter', 'latex')
%         set(gca, 'FontSize', 15)
%         print(['OrderCost-Res01-', num2str(jCase)], '-dpdf')
% 
% 
%         figure(30); 
%         semilogy( ErrorNorms.IterEnd, [col, '*-.'], 'linewidth', 2);
%         hold on
%         xlabel('Iteration', 'interpreter', 'latex')
%         ylabel('NormResidual', 'interpreter', 'latex')
%         set(gca, 'FontSize', 15)
%         print(['OrderCost-ResE-', num2str(jCase)], '-dpdf')
% 
%         figure(31); 
%         loglog( ErrorNorms.IterEnd(1:end-1), ErrorNorms.IterEnd(2:end), [col, '*-.'], 'linewidth', 2);
%         hold on
%         xlabel('$R_i$', 'interpreter', 'latex')
%         ylabel('$R_{i+1}$', 'interpreter', 'latex')
%         set(gca, 'FontSize', 15)
%         print(['OrderCost-ResE1-', num2str(jCase)], '-dpdf')



        iCase = iCase+1;
    end
    jCase = jCase+1;
end

