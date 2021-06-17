function []=main10()

% EFECT OF SIZE
ProblemData = CreateDefaultStructure();

ProblemData.tFinal =  1e-3;
ProblemData.nTimeSteps = 500;
ProblemData.nu = 1/3;



NODES = linspace(1,3,100);
NODES = floor(10.^NODES);
NODES = [4,5,6,7,8,9, NODES];
NODES = unique(NODES);

for i = [1:6]
    figure(i)
    hold off
end

ProblemData.StabMethod=1;
j = 1;

for StabMethod = [1:5]
    i = 1;
    for nodes = NODES
        if ( i > 2)
            ProblemData(j,i) = ProblemData(j,i-1);
        else
            ProblemData(j,i) = ProblemData(1,1);
        end
        ProblemData(j,i).StabMethod = StabMethod;
        ProblemData(j,i).PrimalForm = false;
        ProblemData(j,i).nNodes = nodes;
        ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
        i = i+1;
    end
    SPEC = '';
    if ( StabMethod == 1)
        SPEC = 'k*-';
    elseif ( StabMethod == 3)
        SPEC = '--';
    elseif ( StabMethod == 5)
        SPEC = 'g-.';
    end
    PlotTheDataAndAddLegend('$h_e$ (m)', [1./([ProblemData(j,:).nNodes]-1)], ProblemData(j,:),SPEC);
    j = j+1;
end

for i = 1:6
    figure(i)
    grid minor
    set(gca, 'FontSize', 14)
    yy = ylim();
    PD = ProblemData(1,1);
    dt = PD.tFinal/PD.nTimeSteps;
    M = PD.M;
    k = PD.k;
    beta = 0.3025;
    gamma = 0.6;
    hCrit = sqrt(6*M*beta*dt*k/gamma );
     
    plot( hCrit*[1,1], yy, 'k:')
    
    ll = legend('P1-P1-P1','White \& Borja', 'Sun et al', 'Li \& Wei', 'This work', 'location', 'best');
    set(ll, 'interpreter', 'latex');
     print(['NFigure_', num2str(i)], '-dpdf')
    hold off
end



% EFECT OF TIME
% ProblemData = CreateDefaultStructure();
% ProblemData.PrimalForm = true;
% ProblemData.tFinal =  1e-9;
% ProblemData.nTimeSteps = 10;
% ProblemData.nu = 1/3;
% ProblemData.nNodes = 100;
% 
% M = linspace(1,8,30);
% M = (10.^M);
% 
% for i = [1:6]
%     figure(i)
%     hold off
% end
% j = 1;
% for PrimalMixed = [1,0,1,0]
%     i = 1;
%     for m = M
%         if ( i > 2)
%             ProblemData(j,i) = ProblemData(j,i-1);
%         else
%             ProblemData(j,i) = ProblemData(1,1);
%         end
%         ProblemData(j,i).StabMethod = 1;
%         if (j > 2)
%             ProblemData(j,i).StabMethod=5;
%         end
%         ProblemData(j,i).M = m;
%         ProblemData(j,i).k = 1/m;
%         ProblemData(j,i).PrimalForm = PrimalMixed;
%         ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
%         i = i+1;
%     end
%     PlotTheDataAndAddLegend('$M$ (kPa)', M, ProblemData(j,:));
%     j = j+1;
% end
% 
% for i = 1:4
%     figure(i)
%     yy = ylim();
% %     plot(( 1/((ProblemData(1,1).nNodes-1)^2)/6)*[1,1], yy, 'k:')
% %     print(['EffectTimeSteps_', num2str(i)], '-dpng')
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
ylabel('$\|p(x_i) - p^h_i\|_{\infty}$ (kPa)','interpreter', 'latex')
pause(0.00001)

figure(3)
loglog( xData, [ProblemData(:).L2DISP],  SPEC, 'linewidth', SSS)
ylabel('$\| u-u_h\|_{L2}^*$ (m)', 'interpreter', 'latex')
pause(0.00001)

figure(4)
loglog( xData, [ProblemData(:).LInfDISP],  SPEC, 'linewidth', SSS)
ylabel('$\|u(x_i) - u^h_i\|_{\infty}$ (m)','interpreter', 'latex')
pause(0.0001)

figure(5)
loglog( xData, [ProblemData(:).L2WATER],  SPEC, 'linewidth', SSS)
ylabel('$\| w-w_h\|_{L2}^*$ (m)', 'interpreter', 'latex')
pause(0.00001)

figure(6)
loglog( xData, [ProblemData(:).LInfWATER],  SPEC, 'linewidth', SSS)
ylabel('$\|w(x_i) - w^h_i\|_{\infty}$ (m)','interpreter', 'latex')
pause(0.0001)

for fig = 1:6
    figure(fig)
%     ll = legend('primal','mixed','primal-S', 'mixed-S', 'location', 'best');
    %ll = legend('P1-P1','P2-P1','P3-P2','P4-P3','P1-P1-Stab', 'location','best');
    %ll = legend('P1-P1-P1','White \& Borja', 'Sun et al', 'Li \& Wei', 'This work', 'location', 'best');
    %set(ll, 'interpreter', 'latex');
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
% ProblemData.M = 10^0;
% ProblemData.k = 10^0;
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
ProblemData.LInfWATER = 0;
ProblemData.L2WATER = 0;

ProblemData.TT = 0;
ProblemData.XX = 0;

ProblemData.xxPW = [];
ProblemData.xxU = [];
ProblemData.OrderU = 1;
ProblemData.OrderPW = 1;

