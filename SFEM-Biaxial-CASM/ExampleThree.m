function [] = ExampleThree()
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 20;

figure(212); clf;
figure(214); clf;

for eSize= [0.25, 0.15, 0.10]
    
    if ( eSize == 0.10)
        XNAME = 'Fine';
        SPEC = '-.';
    elseif ( eSize == 0.15)
        XNAME = 'Int';
        SPEC = ':';
    else 
        XNAME = 'Coarse';
        SPEC = '-';
    end
    
    model = createpde(1);
    
    H = -2;
    R1 = [3,4, 0, 1, 1, 0, H, H, 0, 0]';
    
    
    
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
    
    
    
    nSteps = 5000;
    dt = 0.03/nSteps;
    
    
    
    
    tic
    [U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
    toc
    FF = [information2.F];
    FF(1:2:end) = FF(1:2:end);
    figure(212); 
    plot( [information2.t], FF(1:2:end), ['r', SPEC], 'linewidth', 2,'DisplayName', ['NS-T3T3. ' XNAME])
    hold on
    print('Biaxial-Reaction', '-dpdf')


    figure(214); 
    plot( [information2.t], FF(2:2:end), ['r', SPEC] , 'linewidth', 2, 'DisplayName', ['NS-T3T3. ' XNAME])
    hold on
    print('Biaxial-Water', '-dpdf');


    
    figure(557); clf
    pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
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
    
    
    
    figure(457); clf
    GPNodes = ComputeStrainInvatiants(GPNodes);
    PlotHistoryVariableNodal( Nodes1, Elements1, GPNodes, [GPNodes.StrainDev]');
    drawnow
    
    
    
    
    
    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
    toc
    FF = [information.F];
    FF(1:2:end) = FF(1:2:end);
    figure(212)
    plot( [information.t], FF(1:2:end), ['g', SPEC], 'linewidth', 2,'DisplayName',['T3T3. ' XNAME])
    print('Biaxial-Reaction', '-dpdf')
    figure(214)
    plot( [information.t], FF(2:2:end), ['g', SPEC], 'linewidth', 2,'DisplayName', ['T3T3. ' XNAME])
    print('Biaxial-Water', '-dpdf');
    
    figure(556); clf
    pdeplot(model1,'XYData',U(3:3:end),'ColorMap','jet');
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
    
    figure(456); clf
    GPInfo = ComputeStrainInvatiants(GPInfo);
    PlotHistoryVariable( Nodes1, Elements1, GPInfo, [GPInfo.StrainDev]');
    drawnow
    
    
    
    
    
    
    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
    toc
    
    FF = [information.F];
    FF(1:2:end) = FF(1:2:end);
    figure(212)
    plot( [information.t], FF(1:2:end), ['b', SPEC], 'linewidth', 2, 'DisplayName',  ['T6T3. ' XNAME])
    hold on
    print('Biaxial-Reaction', '-dpdf')
    
    figure(214)
    plot( [information.t], FF(2:2:end), ['b', SPEC], 'linewidth', 2, 'DisplayName', ['T6T3. ' XNAME])
    hold on
    print('Biaxial-Water', '-dpdf');
    
    figure(559); clf
    pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet');
    drawnow
    
    figure(959); clf
    SV = [];
    pEff = [];
    StrainDev = [];
    GPInfo = ComputeStrainInvatiants(GPInfo);
    for i = 1:size(GPInfo,1)
        for j = 1:size(GPInfo, 2)
            SV(i,j) = GPInfo(i,j).StressNew(2);
            pEff(i,j) = mean(GPInfo(i,j).StressNew(1:3));
            StrainDev(i,j) = GPInfo(i,j).StrainDev;
        end
    end
    PlotHistoryVariable( Nodes, Elements, GPInfo, SV);
    drawnow
    
    
    figure(359); clf
    PlotHistoryVariable( Nodes, Elements, GPInfo, pEff);
    drawnow
    
    figure(459); clf
    PlotHistoryVariable( Nodes, Elements, GPInfo, StrainDev);
    drawnow
    
    figure(957)
    cc = caxis;
    i = 1;
    pause(1)
    for iii = [956, 957, 959]
        figure(iii)
        axis equal; xlim([0,1]); ylim([H, 0]); axis off
        colormap jet
        caxis(cc);
        colorbar
        drawnow
        pause(1)
        
        fig = figure(iii);
        exportgraphics(fig,['Biaxial-SV-', XNAME, num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    
    figure(357)
    cc = caxis;
    i = 1;
    pause(1)
    for iii = [356, 357, 359]
        figure(iii);
        axis equal; xlim([0,1]); ylim([H, 0]); axis off
        colormap jet
        caxis(cc);
        colorbar
        drawnow
        pause(1)
        fig = figure(iii);
        exportgraphics(fig,['Biaxial-pEff-', XNAME, num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    
    figure(557)
    cc = caxis;
    i = 1;
    pause(1)
    for iii = [556, 557, 559]
        figure(iii)
        axis equal; xlim([0,1]); ylim([H, 0]); axis off
        colormap jet
        caxis(cc);
        colorbar
        drawnow
        pause(1)
        fig = figure(iii);
        exportgraphics(fig,['Biaxial-Water-', XNAME, num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    
    
    figure(457)
    cc = caxis;
    i = 1;
    pause(1)
    for iii = [456, 457, 459]
        figure(iii)
        axis equal; xlim([0,1]); ylim([H, 0]); axis off
        colormap jet
        caxis(cc);
        colorbar
        drawnow
        pause(1)
        fig = figure(iii);
        exportgraphics(fig,['Biaxial-DevStrain-', XNAME, num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    figure(212)
    legend('location', 'best', 'interpreter', 'latex')
    set(gca, 'FontSize', 15)
    xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
    ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
    fig = figure(212);
    exportgraphics(fig,'Biaxial-Reaction.pdf', 'BackgroundColor', 'none','ContentType','vector');
    
    figure(214)
    legend('location', 'best', 'interpreter', 'latex')
    set(gca, 'FontSize', 15)
    xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
    ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
    fig = figure(214);
    exportgraphics(fig,'Biaxial-Water.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
    
end


function GPInfo = ComputeStrainInvatiants(GPInfo)

for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        Strain = GPInfo(i,j).StrainNew;
        eVol = sum(Strain(1:3));
        eDev = 0;
        for k = 1:3
            eDev = eDev + (Strain(k)-eVol/3)^2;
        end
        for k = 4:6
            eDev = eDev + 2*(1/2*Strain(k))^2;
        end
        GPInfo(i,j).StrainVol = eVol;
        GPInfo(i,j).StrainDev = sqrt(eDev);
    end
end


