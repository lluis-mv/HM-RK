function []=main()


% Effect of the size
ProblemData = CreateDefaultStructure();
% close all; clf;

ProblemData.PrimalForm = true;

ProblemData.tFinal =  0.95;
ProblemData.nTimeSteps =  1;
NNODES = linspace(1,3.5,20);
NNODES = 10.^NNODES;
NNODES = floor(NNODES)

NNODES = [2,3,4,5,6,7,8,9,NNODES];
NNODES = linspace(1,2.6,10);
NNODES = 10.^NNODES;
NNODES = floor(NNODES)

NNODES = [2,3,4,5,6,7,8,9,NNODES];
for i = [1,2,3,4,101,102,103, 104]
    figure(i)
    hold off
end


j = 1;
for StabMethod = [1,2,3,4]
    i = 1;
    for nNodes = NNODES
        if ( i > 2)
            ProblemData(j,i) = ProblemData(j,i-1);
        else
            ProblemData(j,i) = ProblemData(1,1);
        end
        ProblemData(j,i).nNodes = nNodes;
        ProblemData(j,i).StabMethod = StabMethod;
        ProblemData(j,i).OrderU = StabMethod;
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
    %plot( sqrt(6*ProblemData(1,1).tFinal/ProblemData(1,1).nTimeSteps)*[1,1], [1e-2,1e-1], ':k')
end






function []=PlotTheDataAndAddLegend(XLABEL, xData, ProblemData,SPEC)
if (nargin == 3)
    SPEC = '';
end
figure(1)
loglog( xData, [ProblemData(:).L2],  ['*-',SPEC])
ylabel('$\| p-p_h\|_{L2}$', 'interpreter', 'latex')
pause(0.0001)

figure(2)
loglog( xData, [ProblemData(:).LInf],  ['*-',SPEC])
ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
pause(0.00001)

figure(3)
loglog( xData, [ProblemData(:).L2DISP],  ['*-',SPEC])
ylabel('$\| u-u_h\|_{L2}$', 'interpreter', 'latex')
pause(0.00001)

figure(4)
loglog( xData, [ProblemData(:).LInfDISP],  ['*-',SPEC])
ylabel('$\|u(x_i) - u^h_i\|$ (kPa)','interpreter', 'latex')
pause(0.0001)

for fig = 1:4
    figure(fig)
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', '$P2-P1$','$P3-P1$','location', 'best');
    ll = legend('P1-P1','P2-P1','P3-P2','P4-P3','P1-P1-Stab', 'location','best');
    set(ll, 'interpreter', 'latex');
    grid minor
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
ProblemData.QBiot = 1e64;
ProblemData.M = 1;
ProblemData.k = 1;
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

