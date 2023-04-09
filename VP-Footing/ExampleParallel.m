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
% MakeSketch(eSize)
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

nSteps = 10;



ind = find(Nodes(:,2) == max( Nodes(:,2)));
xx = sort(Nodes(ind,1));
ind = find(xx == 1);
l = 0.5*(xx(ind)+xx(ind+1));
l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));



NPARA = [1,2,4,8];
TIME = nan*ones( 8, length(NPARA));

color = '';

for RK = [1,2,4,6,8]
    i = 1;
    for Npara = NPARA

        ETA = 100;

        CP.PerzynaEta = ETA;
        CP.MCC = 3;
        CP.RK = -RK;
        if ( CP.PerzynaEta == 0)
            CP.MCC = 4;
        end


        myCluster=parcluster('local'); myCluster.NumWorkers=Npara; parpool(myCluster,Npara)

        tic
        [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
        TIME(RK, i) = toc;
        i = i+1;
        delete(gcp('nocreate'));

        FF = [information.F];
        FF(1:2:end) = FF(1:2:end)/l2;
        figure(212)
        plot( [information.t]*indentation, FF(1:2:end), color, 'linewidth', 2, 'DisplayName',  ['$\eta = $', num2str(ETA)])
        hold on


        figure(214)
        plot( [information.t]*indentation, FF(2:2:end), color, 'linewidth', 2, 'DisplayName',  ['$\eta = $', num2str(ETA)])
        hold on

        figure(21); clf;
        plot(NPARA, TIME./TIME(1,1), '*-.')
        
        xlabel('Number of cores', 'interpreter', 'latex')
        ylabel('Computational cost (s)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        drawnow

        print('Parallel1', '-dpdf')

        figure(21); clf;

        for jj = 1:size(TIME,1)
            TIME2(jj,:) = TIME(jj,:)./TIME(jj,1);
        end
        plot(NPARA, TIME2, '*-.')
        
        xlabel('Number cores', 'interpreter', 'latex')
        ylabel('Computational cost (s)', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        drawnow

        print('Parallel2', '-dpdf')

    end
end
