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


model = createpde(1);


R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';



ESIZE = 1./[2:9];



nSteps = 100;
dt = 1.0/nSteps;


i = 1;

g = decsg(R1);
geometryFromEdges(model, g);

model1 = createpde(1);
geometryFromEdges(model1, g);


if (true)

    for eSize = ESIZE

        mesh = generateMesh(model, 'Hmax', eSize);
        Nodes = mesh.Nodes';
        Elements = mesh.Elements';



        mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodes1 = mesh1.Nodes';
        Elements1 = mesh1.Elements';

        figure(1)
        triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')
        axis equal;
        axis off;
        print(['Mesh-', num2str(i)], '-dpdf');

        [GPInfo] = ComputeElementalMatrices(Nodes1, Elements1, CP, 'T3T3');
        he = mean(sqrt([GPInfo(:,:).Weight]));

        eSizeAxis(i) = he;



        tic
        [U, GPInfo, GPNodes, rrr,  information2, nZero] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
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
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
        TIMElinear(i)= toc;
        nZerolinear(i) = nZero;
        TIMEnodal./TIMElinear
        FF = [information.F];
        PWlinear(i) = FF(end);
        Qlinear(i) = FF(end-1)/l;

        tic
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', 1);
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
            'TIMEquad', 'nDofsquad', 'PWquad', 'Qquad', 'nZeroquad');


        figure(99); clf
        plot(eSizeAxis, Qnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(eSizeAxis, Qlinear, 'gs-.', 'DisplayName', 'T3T3')
        plot(eSizeAxis, Qquad, 'bs-.', 'DisplayName', 'T6T3')
        drawnow
        xlabel('$h_e$ (m)', 'interpreter', 'latex')
        ylabel('Footing resistance (kPa)', 'interpreter', 'latex')

        figure(100); clf
        plot(eSizeAxis, PWnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(eSizeAxis, PWlinear, 'gs-.', 'DisplayName', 'T3T3')
        plot(eSizeAxis, PWquad, 'bs-.', 'DisplayName', 'T6T3')
        xlabel('$h_e$ (m)', 'interpreter', 'latex')
        ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
        drawnow

        figure(101); clf
        plot(nDofs, TIMEnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(nDofs, TIMElinear, 'g*-.', 'DisplayName', 'T3T3')
        plot(nDofsquad, TIMEquad, 'b*-.', 'DisplayName', 'T6T3')
        drawnow
        xlabel('Number of dofs', 'interpreter', 'latex')
        ylabel('Computational cost (s)', 'interpreter', 'latex')
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')


        figure(102); clf
        plot(nDofs, TIMEnodal./nDofs, 'r*-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(nDofs, TIMElinear./nDofs, 'g*-.', 'DisplayName', 'T3T3')
        plot(nDofsquad, TIMEquad./nDofsquad, 'b*-.', 'DisplayName', 'T6T3')
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

        figure(99); drawnow; pause(1); print('Footing-Load-1', '-dpdf')
        figure(100); drawnow; pause(1); print('Footing-WP-1', '-dpdf')
        figure(101); drawnow; pause(1); print('Footing-Cost-1', '-dpdf')
        figure(102); drawnow; pause(1); print('Footing-Velocity-1', '-dpdf')
        i = i+1;
    end

end



eSize = 1/4;
mesh = generateMesh(model, 'Hmax', eSize);
Nodes = mesh.Nodes';
Elements = mesh.Elements';



mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
Nodes1 = mesh1.Nodes';
Elements1 = mesh1.Elements';

NSTEPS = 10.^linspace(1.5,log10(1000),10); NSTEPS(end) = 1000;
NSTEPS = floor(NSTEPS);
i = 1;


for nSteps = NSTEPS


    dt = 1.0/nSteps;
    dtAxis(i) = dt;



    tic
    [U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);

    TIMEnodal(i)= toc;
    nDofs(i) = size(Nodes1,1)*3;
    ind = find(Nodes(:,2) == max( Nodes(:,2)));
    xx = sort(Nodes(ind,1));
    ind = find(xx == 1);
    l = 0.5*(xx(ind)+xx(ind+1));
    l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));

    FF = [information2.F];
    PWnodal(i) = FF(end);
    Qnodal(i) = FF(end-1)/l;



    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
    TIMElinear(i)= toc;
    TIMEnodal./TIMElinear
    FF = [information.F];
    PWlinear(i) = FF(end);
    Qlinear(i) = FF(end-1)/l;

    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', 1);
    TIMEquad(i) = toc;
    FF = [information.F];
    PWquad(i) = FF(end);
    Qquad(i) = FF(end-1)/l2;
    nDofsquad(i) = size(Nodes,1)*3;



    figure(99); clf
    plot(NSTEPS(1:i), Qnodal(1:i), 'rs-.', 'DisplayName', 'NS-T3T3')
    hold on
    plot(NSTEPS(1:i), Qlinear(1:i), 'gs-.', 'DisplayName', 'T3T3')
    plot(NSTEPS(1:i), Qquad(1:i), 'bs-.', 'DisplayName', 'T6T3')
    drawnow
    xlabel('$n_{steps}$', 'interpreter', 'latex')
    ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
    set(gca, 'XScale', 'log')

    figure(100); clf
    plot(NSTEPS(1:i), PWnodal(1:i), 'rs-.', 'DisplayName', 'NS-T3T3')
    hold on
    plot(NSTEPS(1:i), PWlinear(1:i), 'gs-.', 'DisplayName', 'T3T3')
    plot(NSTEPS(1:i), PWquad(1:i), 'bs-.', 'DisplayName', 'T6T3')
    xlabel('$n_{steps}$', 'interpreter', 'latex')
    ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
    drawnow
    set(gca, 'XScale', 'log')

    save('TimeUndrainedData.mat', ...
        'i', 'dtAxis', 'NSTEPS', ...
        'PWnodal', 'Qnodal',  ...
        'PWlinear', 'Qlinear', ...
        'PWquad', 'Qquad');

    for jj = [99, 100]
        figure(jj)
        set(gca, 'FontSize', 15)
        ll = legend();
        set(ll, 'location', 'best', 'interpreter', 'latex')
    end

    figure(99); drawnow; pause(1); print('Footing-Load-2', '-dpdf')
    figure(100); drawnow; pause(1); print('Footing-WP-2', '-dpdf')
    i = i+1;
end
