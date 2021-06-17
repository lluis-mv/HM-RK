function []=main6()
clf; close all; clear all; clc;
if ( false)
    ProblemData = CreateDefaultStructure();
    
    ProblemData.tFinal =  1e-5;
    ProblemData.nTimeSteps =  1;
    ProblemData.PrimalForm = true;
    
    
    
    NU = linspace(-4,0,20)
    NU = 0.5*(1-10.^NU)
    j = 1;
    for StabMethod = [1:5]
        i = 1;
        for nu = NU
            if ( i > 2)
                ProblemData(j,i) = ProblemData(j,i-1);
            else
                ProblemData(j,i) = ProblemData(1,1);
            end
            ProblemData(j,i).nu = nu;
            ProblemData(j,i).StabMethod = StabMethod;
            
            ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
            i = i+1;
        end
        PlotTheDataAndAddLegend('0.5-$\nu$', 0.5-[ProblemData(j,:).nu], ProblemData(j,:));
        j = j+1;
    end
    
    
    % The same but with a much higher dt
    ProblemData = CreateDefaultStructure();
    close all; clf;
    ProblemData.tFinal =  1e-4;
    ProblemData.nTimeSteps =  1;
    ProblemData.PrimalForm = true;
    
    
    
    NU = linspace(-4,0,20)
    NU = 0.5*(1-10.^NU)
    j = 1;
    for StabMethod = [1:5]
        i = 1;
        for nu = NU
            if ( i > 2)
                ProblemData(j,i) = ProblemData(j,i-1);
            else
                ProblemData(j,i) = ProblemData(1,1);
            end
            ProblemData(j,i).nu = nu;
            ProblemData(j,i).StabMethod = StabMethod;
            
            ProblemData(j,i) = SolveProblemAndComputeNorms( ProblemData(j,i));
            i = i+1;
        end
        PlotTheDataAndAddLegend('0.5-$\nu$', 0.5-[ProblemData(j,:).nu], ProblemData(j,:));
        j = j+1;
    end
    
end

% Effect of the size
ProblemData = CreateDefaultStructure();

ProblemData.tFinal =  1e-5;
ProblemData.nTimeSteps =  1;
ProblemData.PrimalForm = true;

ProblemData.tFinal =  1e-9;
NNODES = linspace(1,3,20);
NNODES = 10.^NNODES;
NNODES = floor(NNODES)
j = 1;
for StabMethod = [1,6,7,8]
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
    PlotTheDataAndAddLegend('$h_e$ (m)', 1./([ProblemData(j,:).nNodes]-1), ProblemData(j,:));
    j = j+1;
end

for i = 1:4
    figure(i)
    plot( sqrt(6*ProblemData(1,1).tFinal)*[1,1], [1e-2,1e-1], ':k')
end






function []=PlotTheDataAndAddLegend(XLABEL, xData, ProblemData)

figure(1)
loglog( xData, [ProblemData(:).L2],  '*-')
ylabel('$\| p-p_h\|_{L2}$', 'interpreter', 'latex')
pause(0.0001)

figure(2)
loglog( xData, [ProblemData(:).LInf],  '*-')
ylabel('$\|p(x_i) - p^h_i\|$ (kPa)','interpreter', 'latex')
pause(0.00001)

figure(3)
loglog( xData, [ProblemData(:).L2DISP],  '*-')
ylabel('$\| u-u_h\|_{L2}$', 'interpreter', 'latex')
pause(0.00001)

figure(4)
loglog( xData, [ProblemData(:).LInfDISP],  '*-')
ylabel('$\|u(x_i) - u^h_i\|$ (kPa)','interpreter', 'latex')
pause(0.0001)

for fig = 1:4
    figure(fig)
    ll = legend('Non-stabilized', 'White \& Borja ($1/(2G)$)', 'Sun et al', 'Li \& Wei ($3/M$)', 'This work ($2/M - 12 \Delta t k/ h^2 $)', '$P2-P1$','$P3-P1$','location', 'best');
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