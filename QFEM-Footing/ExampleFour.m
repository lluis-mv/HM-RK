function [] = ExampleFour()
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
CP.MCC = 2;


ESIZE = [2:11];



nSteps = 100;
dt = 1.0/nSteps;



i = 1;



if (true)

    for eSize = ESIZE


        [Nodes, Elements] = ReadTheMesh(['Mesh', num2str(eSize), '.msh']);
        [Nodes1, Elements1] = SimplifyOrder(Nodes,Elements);

        

        fig = figure(1)
        PlotTheMesh( Elements1, Nodes1)
        %triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')
        axis equal;
        axis off;
        exportgraphics(fig,['MeshQ-', num2str(i), '.pdf'], 'BackgroundColor', 'none','ContentType','vector');

        [GPInfo] = ComputeElementalMatrices(Nodes1, Elements1, CP, 'Q4Q4');
        he = mean(sqrt([GPInfo(:,:).Weight]));

        eSizeAxis(i) = he;



        tic
        [U, GPInfo, GPNodes, rrr,  information2, nZero] = ComputeImplicitNonLinearProblemNodalQuad(Nodes1, Elements1, CP, dt, nSteps, 'Q4Q4', 1);
        TIMEnodal(i)= toc;
        nDofs(i) = size(Nodes1,1)*3;
        nZeronodal(i) = nZero;
        ind = find(Nodes(:,2) == max( Nodes(:,2)));
        xx = sort(Nodes(ind,1));
        ind = find(xx == 1);
        l = 0.5*(xx(ind)+xx(ind+1));
        l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));


        FF = [information2.F];
        PWnodal(i) = FF(end);
        Qnodal(i) = FF(end-1)/l;


        tic
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'Q4Q4', 1);
        TIMElinear(i)= toc;
        nZerolinear(i) = nZero;
        TIMEnodal./TIMElinear
        FF = [information.F];
        PWlinear(i) = FF(end);
        Qlinear(i) = FF(end-1)/l;


%         tic
%         [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'M3T3', 1);
        TIMEmixed(i)= nan;
        nDofsmixed(i) = size(Nodes1,1)*4*nan;
        nZeromixed(i) = nZero*nan;
        FF = [information.F];
        PWmixed(i) = FF(end)*nan;
        Qmixed(i) = FF(end-1)/l*nan;


        tic
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'Q8Q4', 1);
        TIMEquad(i) = toc;
        nZeroquad(i) = nZero;
        FF = [information.F];
        PWquad(i) = FF(end);
        Qquad(i) = FF(end-1)/l2;
        nDofsquad(i) = size(Nodes,1)*3;


        save('UndrainedData.mat', ...
            'ESIZE', 'i', 'eSizeAxis', ...
            'TIMEnodal', 'nDofs', 'PWnodal', 'Qnodal', 'nZeronodal', ...
            'TIMElinear',  'PWlinear', 'Qlinear', 'nZerolinear', ...
            'TIMEquad', 'nDofsquad', 'PWquad', 'Qquad', 'nZeroquad', ...
            'TIMEmixed', 'nDofsmixed', 'PWmixed', 'Qmixed', 'nZeromixed');


        figure(99); clf
        plot(eSizeAxis, Qnodal, 'rs-.', 'DisplayName', 'NS-Q4Q4')
        hold on
        plot(eSizeAxis, Qlinear, 'gs-.', 'DisplayName', 'Q4Q4')
        
        plot(eSizeAxis, Qquad, 'bs-.', 'DisplayName', 'Q8Q4')
        drawnow
        xlabel('$h_e$ (m)', 'interpreter', 'latex')
        ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
        ylim([32,37])
        drawnow

        figure(100); clf
        plot(eSizeAxis, PWnodal, 'rs-.', 'DisplayName', 'NS-Q4Q4')
        hold on
        plot(eSizeAxis, PWlinear, 'gs-.', 'DisplayName', 'Q4Q4')
        
        plot(eSizeAxis, PWquad, 'bs-.', 'DisplayName', 'Q8Q4')
        xlabel('$h_e$ (m)', 'interpreter', 'latex')
        ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
        drawnow
        ylim([18,23])

        drawnow
        figure(101); clf
        plot(nDofs, TIMEnodal, 'r*-.', 'DisplayName', 'NS-Q4Q4')
        hold on
        plot(nDofs, TIMElinear, 'g*-.', 'DisplayName', 'Q4Q4')
        
        plot(nDofsquad, TIMEquad, 'b*-.', 'DisplayName', 'Q8Q4')
        drawnow
        xlabel('Number of dofs', 'interpreter', 'latex')
        ylabel('Computational cost (s)', 'interpreter', 'latex')
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')


        figure(102); clf
        plot(nDofs, TIMEnodal./nDofs, 'r*-.', 'DisplayName', 'NS-Q4Q4')
        hold on
        plot(nDofs, TIMElinear./nDofs, 'g*-.', 'DisplayName', 'Q4Q4')
        
        plot(nDofsquad, TIMEquad./nDofsquad, 'b*-.', 'DisplayName', 'Q8Q4')
        drawnow
        xlabel('Number of dofs', 'interpreter', 'latex')
        ylabel('Computational cost / Number of dofs (s)', 'interpreter', 'latex')
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')

        for jj = [99, 100, 101, 102]
            figure(jj)
            set(gca, 'FontSize', 15)
            ll = legend();
            set(ll, 'location', 'best', 'interpreter', 'latex')
        end

        figure(99); drawnow; pause(1); print('FootingQ-Load-1', '-dpdf')
        figure(100); drawnow; pause(1); print('FootingQ-WP-1', '-dpdf')
        figure(101); drawnow; pause(1); print('FootingQ-Cost-1', '-dpdf')
        figure(102); drawnow; pause(1); print('FootingQ-Velocity-1', '-dpdf')
        i = i+1;
        close all;
        ExampleFour2()
        close all;
    end

end


function [] = PlotTheMesh( Elements, Nodes)

ii = [1,2,3,4,1];
for elem = 1:size(Elements,1)
    Celem = Elements(elem,:);
    Xe = Nodes(Celem,:);
    plot(Xe(ii,1), Xe(ii,2), 'k')
    hold on
end

