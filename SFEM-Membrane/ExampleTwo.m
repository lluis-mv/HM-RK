function [] = ExampleTwo()
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.0;
CP.k = 1E-16;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.Elastic = true;


model = createpde(1);


R1 = [3,4,0, 48, 48, 0, 0, 44, 44+16, 44]';
R1 = [3,4,0, 1, 1, 0, 0, 0, 1, 1]';



ESIZE = 1./[5:3:26];



nSteps = 1;
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
        Nodes = CreateMapping(Nodes);
        
        
        mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodes1 = mesh1.Nodes';
        Elements1 = mesh1.Elements';
        Nodes1 = CreateMapping(Nodes1);

        figure(1)
        triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')
        
        axis equal;
        axis off;
        print(['Mesh-Membrane-', num2str(i)], '-dpdf');
        
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
        
        
        
        
        FF = [information2.F];
        PWnodal(i) = FF(end);
        Qnodal(i) = FF(end-1);
        
        
        tic
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
        TIMElinear(i)= toc;
        nZerolinear(i) = nZero;
        TIMEnodal./TIMElinear
        FF = [information.F];
        PWlinear(i) = FF(end);
        Qlinear(i) = FF(end-1);
        
        tic
        [U, GPInfo, rrr,  information, nZero] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', 1);
        TIMEquad(i) = toc;
        nZeroquad(i) = nZero;
        FF = [information.F];
        PWquad(i) = FF(end);
        Qquad(i) = FF(end-1);
        nDofsquad(i) = size(Nodes,1)*3;
        
        
        save('UndrainedData.mat', ...
            'ESIZE', 'i', 'eSizeAxis', ...
            'TIMEnodal', 'nDofs', 'PWnodal', 'Qnodal', 'nZeronodal', ...
            'TIMElinear',  'PWlinear', 'Qlinear', 'nZerolinear', ...
            'TIMEquad', 'nDofsquad', 'PWquad', 'Qquad', 'nZeroquad');
        
        
        figure(99); clf
        plot(eSizeAxis, Qnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(eSizeAxis, Qlinear, 'g*-.', 'DisplayName', 'T3T3')
        plot(eSizeAxis, Qquad, 'b*-.', 'DisplayName', 'T6T3')
        drawnow
        xlabel('$h_e$ (m)', 'interpreter', 'latex')
        ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
        
        figure(100); clf
        plot(eSizeAxis, PWnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
        hold on
        plot(eSizeAxis, PWlinear, 'g*-.', 'DisplayName', 'T3T3')
        plot(eSizeAxis, PWquad, 'b*-.', 'DisplayName', 'T6T3')
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
        
        figure(99);  drawnow; pause(1); print('Membrane-Load-1', '-dpdf')
        figure(100); drawnow; pause(1); print('Membrane-WP-1', '-dpdf')
        figure(101); drawnow; pause(1); print('Membrane-Cost-1', '-dpdf')
        figure(102); drawnow; pause(1); print('Membrane-Velocity-1', '-dpdf')
        i = i+1;
    end
    
end




function Xn = CreateMapping(X)

Xn = 0*X;

Xn(:,1) = 48*X(:,1);

for i = 1:size(X,1)
    x = X(i,1);
    y = X(i,2);
    c0 = 44*x;
    slope = 44*(1-x)+16*x;
    y = c0+slope*y;
    Xn(i,2) = y;
end

