function []=main6()
% clf; close all; clear all; clc;
% First Figure. Effect of Poisson

ProblemData = CreateDefaultStructure();
ProblemData.PrimalForm = true;
ProblemData.tFinal =  1e-4;
ProblemData.nTimeSteps =  1;
ProblemData.nNodes = 20;
NODES = linspace(-5,0,100)
NODES = 0.5*(1-10.^NODES);

for i = [1,2,3,4]
    figure(i)
    hold off
end
j = 1;
for StabMethod = [1, 10, 2:5]
    i = 1;
    for nu = NODES
        if ( i > 2)
            ProblemData(j,i) = ProblemData(j,i-1);
        else
            ProblemData(j,i) = ProblemData(1,1);
        end
        ProblemData(j,i).StabMethod = StabMethod;
        ProblemData(j,i).nu = nu;
        ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
        i = i+1;
    end
    PlotTheDataAndAddLegend('$0.5-\nu$', (0.5-[ProblemData(j,:).nu]), ProblemData(j,:),'');
    j = j+1;
end

for i = 1:4
    figure(i)
    grid minor
    print(['EffectPoisson_', num2str(i)], '-dpng')
    hold off
end


% EFECT OF SIZE
% ProblemData = CreateDefaultStructure();
% ProblemData.PrimalForm = true;
% ProblemData.tFinal =  4e-4;
% ProblemData.nTimeSteps =  10;
% ProblemData.nu = 1/3;
% NODES = linspace(1,3,40);
% NODES = floor(10.^NODES);
% NODES = [4,6,8, NODES];
% NODES = unique(NODES)
% for i = [1,2,3,4]
%     figure(i)
%     hold off
% end
% j = 1;
% for StabMethod = [1, 10, 2:5]
%     i = 1;
%     for nodes = NODES
%         if ( i > 2)
%             ProblemData(j,i) = ProblemData(j,i-1);
%         else
%             ProblemData(j,i) = ProblemData(1,1);
%         end
%         ProblemData(j,i).StabMethod = StabMethod;
%         ProblemData(j,i).nNodes = nodes;
%         ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
%         i = i+1;
%     end
%     PlotTheDataAndAddLegend('$h_e$ (m)', [1./([ProblemData(j,:).nNodes]-1)], ProblemData(j,:),'');
%     j = j+1;
% end
% 
% for i = 1:4
%     figure(i)
%     grid minor
%     yy = ylim();
%     plot( sqrt(6*ProblemData(1,1).tFinal/ProblemData(1,1).nTimeSteps)*[1,1], yy, 'k:')
%     print(['EffectSize_', num2str(i)], '-dpng')
%     hold off
% end



% EFECT OF TIME
% ProblemData = CreateDefaultStructure();
% ProblemData.PrimalForm = true;
% ProblemData.tFinal =  1e-2;
% 
% ProblemData.nu = 1/3;
% ProblemData.nNodes = 40;
% 
% NSTEPS = linspace(0,3,40);
% NSTEPS = floor(10.^NSTEPS);
% NSTEPS = unique(NSTEPS);
% for i = [1,2,3,4]
%     figure(i)
%     hold off
% end
% j = 1;
% for StabMethod = [1, 10, 2:5]
%     i = 1;
%     for nSteps = NSTEPS
%         if ( i > 2)
%             ProblemData(j,i) = ProblemData(j,i-1);
%         else
%             ProblemData(j,i) = ProblemData(1,1);
%         end
%         ProblemData(j,i).StabMethod = StabMethod;
%         ProblemData(j,i).nTimeSteps =  nSteps;
%         ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
%         i = i+1;
%     end
%     PlotTheDataAndAddLegend('$\Delta t$ (s)', [ProblemData(1,1).tFinal./([ProblemData(j,:).nTimeSteps])], ProblemData(j,:));
%     j = j+1;
% end
% 
% for i = 1:4
%     figure(i)
%     yy = ylim();
%     plot(( 1/((ProblemData(1,1).nNodes-1)^2)/6)*[1,1], yy, 'k:')
%     print(['EffectTimeSteps_', num2str(i)], '-dpng')
%     hold off
% end


% EFECT OF TIME (INCORRECT FIGURE
% ProblemData = CreateDefaultStructure();
% ProblemData.PrimalForm = true;
% 
% 
% ProblemData.nu = 1/3;
% ProblemData.nNodes = 40;
% ProblemData.nTimeSteps = 1.0;
% TFINAL = linspace(-10,0,100);
% TFINAL = (10.^TFINAL);
% 
% for i = [1,2,3,4]
%     figure(i)
%     hold off
% end
% j = 1;
% for StabMethod = [1, 10, 2:5]
%     i = 1;
%     for tFinal = TFINAL
%         if ( i > 2)
%             ProblemData(j,i) = ProblemData(j,i-1);
%         else
%             ProblemData(j,i) = ProblemData(1,1);
%         end
%         ProblemData(j,i).StabMethod = StabMethod;
%         ProblemData(j,i).tFinal =  tFinal;
%         ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
%         i = i+1;
%     end
%     PlotTheDataAndAddLegend('$t$ (s)', [ProblemData(j,:).tFinal], ProblemData(j,:), '-');
%     j = j+1;
% end
% 
% for i = 1:4
%     figure(i)
%     grid minor
%     yy = ylim();
%     xlim([TFINAL(1), TFINAL(end)])
%     plot(( 1/((ProblemData(1,1).nNodes-1)^2)/6)*[1,1].*ProblemData(1,1).nTimeSteps, yy, 'k:')
%     print(['EffectTime1_', num2str(i)], '-dpng')
%     hold off
% end
% 


return;


% Effect of the size
ProblemData = CreateDefaultStructure();
% close all; clf;

ProblemData.PrimalForm = true;

ProblemData.tFinal =  1e-2;
ProblemData.nTimeSteps =  1;
NNODES = linspace(1,3.3,20);
NNODES = 10.^NNODES;
NNODES = floor(NNODES)



NNODES = [2,3,4,5,6,7,8,9,NNODES];
for i = [1,2,3,4]
    figure(i)
    hold off
end


j = 1;
for StabMethod = [1,10,11,12,5]
    i = 1;
    for nNodes = NNODES
        if ( i > 2)
            ProblemData(j,i) = ProblemData(j,i-1);
        else
            ProblemData(j,i) = ProblemData(1,1);
        end
        ProblemData(j,i).nNodes = nNodes;
        ProblemData(j,i).StabMethod = StabMethod;
        
        ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
        i = i+1;
    end
    SPEC = [];
    if (StabMethod == 5)
        SPEC = 'k'
    end
    PlotTheDataAndAddLegend('$h_e$ (m)', 1./([ProblemData(j,:).nNodes]-1), ProblemData(j,:), SPEC);
%     xx = [];
%     for abc = 1:size(ProblemData(j,:),2)
%         xx(abc) = length([ProblemData(j,abc).NodalWaterPressure])+length([ProblemData(j,abc).NodalDisplacement])
%     end
%     PlotTheDataAndAddLegend('$nDOFS$', xx, ProblemData(j,:), SPEC);
    j = j+1;
end

for i = 1:4
    figure(i)
end






function []=PlotTheDataAndAddLegend(XLABEL, xData, ProblemData,SPEC)
if (nargin == 3)
    SPEC = '*-';
end
SSS = 1.4;
figure(1)
loglog( xData, [ProblemData(:).L2], SPEC, 'linewidth', SSS)
ylabel('$\| p-p_h\|_{L2}^*$ (kPa)', 'interpreter', 'latex')
pause(0.0001)

figure(2)
loglog( xData, [ProblemData(:).LInf], SPEC, 'linewidth', SSS)
ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
pause(0.00001)

figure(3)
loglog( xData, [ProblemData(:).L2DISP],  SPEC, 'linewidth', SSS)
ylabel('$\| u-u_h\|_{L2}^*$ (m)', 'interpreter', 'latex')
pause(0.00001)

figure(4)
loglog( xData, [ProblemData(:).LInfDISP],  SPEC, 'linewidth', SSS)
ylabel('$\|u(x_i) - u^h_i\|$ (m)','interpreter', 'latex')
pause(0.0001)

for fig = 1:4
    figure(fig)
    ll = legend('P1-P1','P2-P1','White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work', 'location', 'EastOutside');
    %ll = legend('P1-P1','P2-P1','P3-P2','P4-P3','P1-P1-Stab', 'location','best');
    set(ll, 'interpreter', 'latex');
    xlabel(XLABEL, 'interpreter', 'latex')
    hold on
    pause(0.0001)
end

function ProblemData = CreateDefaultStructure()

ProblemData.nNodes = 50;
ProblemData.tFinal = 0.0001;
ProblemData.nTimeSteps = 1;
ProblemData.StabMethod = 1;
ProblemData.H = 1;
ProblemData.PrimalForm = true;
ProblemData.QBiot = 1e444;
ProblemData.M = 10^2;
ProblemData.k = 10^-2;
ProblemData.nu = 1/3;
ProblemData.LOAD = -1;
ProblemData.FixedWP = 0;
ProblemData.NodalWaterPressure = [];
ProblemData.NodalDisplacement = [];
ProblemData.NodalWaterDisplacement=[];
ProblemData.LInfDISP = 0;
ProblemData.L2DISP = 0;
ProblemData.LInf = 0;
ProblemData.L2 = 0;
ProblemData.TT = 0;
ProblemData.XX = 0;

ProblemData.xxPW = [];
ProblemData.xxU = [];
ProblemData.OrderU = 1;
ProblemData.OrderPW = 1;

