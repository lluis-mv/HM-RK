function [] = ExampleThree()
addpath('../')
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

eSize= 1/3;

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



figure(1); clf;
[Elements1, Nodes1] = CreateMesh(0,4, -4, 0, 0.25, 0.25);

triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'b')
[Nodes2] = AntiLaplacianSmoothing(Elements1, Nodes1);



nSteps = 100;
dt = 0.15/nSteps;


for i = 1:6
    figure(211+i); clf;
end



for this = [-0.25, -0.15, 0, 0.25]
    
    Nodes1 = Nodes2;
    
    for i = 1:10
        [Elements1, Nodes1] = AntiLaplacianSmoothing(Elements1, Nodes1, this);
        figure(2); clf;
        triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'r')
    end
    
    figure(21); clf;
    triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')
    axis equal;
    axis off;
    drawnow
    
    
    tic
    [U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
    toc
    ind = find(Nodes(:,2) == max( Nodes(:,2)));
    xx = sort(Nodes(ind,1));
    ind = find(xx == 1);
    l = 0.5*(xx(ind)+xx(ind+1));
    l2 = xx(ind)+0.25*(xx(ind+1)-xx(ind));
    FF = [information2.F];
    FF(1:2:end) = FF(1:2:end)/l;
    figure(212)
    plot( [information2.t], FF(1:2:end), '', 'linewidth', 2,'DisplayName', ['NS-T3T3'])
    hold on
    drawnow
    figure(215)
    plot( [information2.t], FF(2:2:end), '', 'linewidth', 2, 'DisplayName', ['NS-T3T3'])
    hold on
    drawnow
    
    
    figure(557); clf
    pdeplot(Nodes1', Elements1', 'XYData', U(3:3:end), 'ColorMap', 'jet')
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
    
    
    
    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', 1);
    toc
    FF = [information.F];
    FF(1:2:end) = FF(1:2:end)/l;
    figure(213)
    plot( [information.t], FF(1:2:end), '', 'linewidth', 2,'DisplayName', ['T3T3'])
    hold on
    drawnow
    figure(216)
    plot( [information.t], FF(2:2:end), '', 'linewidth', 2,'DisplayName', ['T3T3'])
    hold on
    drawnow
    
    
    figure(556); clf
    pdeplot(Nodes1', Elements1', 'XYData', U(3:3:end), 'ColorMap', 'jet')
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
    
    
    
    [NodesX, ElementsX] = ConvertQuadratic(Nodes1, Elements1);
    
    tic
    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(NodesX, ElementsX, CP, dt, nSteps, 'T6T3');
    toc
    FF = [information.F];
    FF(1:2:end) = FF(1:2:end)/l2;
    figure(214)
    plot( [information.t], FF(1:2:end), '', 'linewidth', 2,'DisplayName', ['T3T3'])
    hold on
    drawnow
    figure(217)
    plot( [information.t], FF(2:2:end), '', 'linewidth', 2,'DisplayName', ['T3T3'])
    hold on
    drawnow
    
    figure(559); clf
    pdeplot(NodesX', ElementsX','XYData',U(3:3:end),'ColorMap','jet');
    drawnow
    
    figure(959); clf
    SV = [];
    pEff = [];
    for i = 1:size(GPInfo,1)
        for j = 1:size(GPInfo, 2)
            SV(i,j) = GPInfo(i,j).StressNew(2);
            pEff(i,j) = mean(GPInfo(i,j).StressNew(1:3));
        end
    end
    PlotHistoryVariable( NodesX, ElementsX, GPInfo, SV);
    drawnow
    
    
    figure(359); clf
    PlotHistoryVariable( NodesX, ElementsX, GPInfo, pEff);
    drawnow
    
    
    
    
    figure(957)
    cc = caxis;
    i = 1;
    pause(1)
    
    thisN = num2str(this);
    index = find(thisN == '.')
    thisN(index) = '_';
    for iii = [956, 957, 959]
        figure(iii)
        axis equal; xlim([0,4]); ylim([-4, 0]); axis off
        colormap jet
        caxis(cc);
        colorbar
        drawnow
        pause(1)
        
        fig = figure(iii);
        exportgraphics(fig,['F2-SV-', num2str(i), '-', thisN,'.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        
        i = i+1;
    end
    
    
    figure(357)
    i = 1;
    pause(1)
    for iii = [356, 357, 359]
        figure(iii)
        axis equal; xlim([0,4]); ylim([-4, 0]); axis off
        colormap jet
        caxis([-13,-5]);
        colorbar
        drawnow
        pause(1)
       
        fig = figure(iii);
        exportgraphics(fig,['F2-pEff-', num2str(i), '-', thisN, '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    
    figure(557)
    i = 1;
    pause(1)
    for iii = [556, 557, 559]
        figure(iii)
        axis equal; xlim([0,4]); ylim([-4, 0]); axis off
        colormap jet
        caxis([0,19]);
        colorbar
        drawnow
        pause(1)
        print(['F2-Water-', num2str(i), '-', thisN], '-dpdf');
        fig = figure(iii);
        exportgraphics(fig,['F2-SV-', num2str(i), '-', thisN, '.pdf'], 'BackgroundColor', 'none','ContentType','vector');
        i = i+1;
    end
    
    
    figure(21)
    pause(1)
    print(['F2-Mesh-', thisN], '-dpdf');
    
    
    for i = 1:3
        fig = figure(211+i)
        ll = legend('Mesh A', 'Mesh B', 'Mesh C', 'Mesh D');
        set(ll, 'location', 'best', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
        ylabel('Footing reaction (kPa)', 'interpreter', 'latex')
        exportgraphics(fig, ['F2-Reaction-', num2str(i), '.dpdf'], 'BackgroundColor', 'none','ContentType','vector');
        
        
        fig = figure(214+i)
        ll = legend('Mesh A', 'Mesh B', 'Mesh C', 'Mesh D');
        set(ll, 'location', 'best', 'interpreter', 'latex')
        set(gca, 'FontSize', 15)
        xlabel('Footing indentation, $u_z/R$', 'interpreter', 'latex')
        ylabel('Water pressure, $p_w$ (kPa)', 'interpreter', 'latex')
        exportgraphics(fig, ['F2-Water-', num2str(i), '.dpdf'], 'BackgroundColor', 'none','ContentType','vector');
    end
end

function [C, X ] = AntiLaplacianSmoothing(C, X, this)
nNodes = size(X,1);
IsNeig = zeros(nNodes,1);

ind = find(X(:,1) == min(X(:,1)));
IsNeig(ind) = 1;

ind = find(X(:,1) == max(X(:,1)));
IsNeig(ind) = 1;
ind = find(X(:,2) == min(X(:,2)));
IsNeig(ind) = 1;
ind = find(X(:,2) == max(X(:,2)));
IsNeig(ind) = 1;
xL = X;
if (nargout == 1)
    for i = 1:nNodes
        if ( IsNeig(i) ==0)
            X(i,:) =  X(i,:) + 0.05*(rand(1,2)-0.5);
        end
    end
    C = X;
    return;
end

for i = 1:nNodes
    if ( IsNeig(i) ==0)
        [a,b] = find(C == i);
        %         NeigNodes = unique(C(a,:));
        NeigNodes = sort(C(a,:));
        [index] = find( NeigNodes~=i);
        NeigNodes = NeigNodes(index);
        xL(i,:) = mean(X(NeigNodes,:));
    end
end
displ = xL -X;
for i = 1:nNodes
    if ( norm(displ(i,:) )> 0)
        X(i,:) = X(i,:) +this *displ(i,:); %/norm(displ(i,:));
    end
end
C = delaunay(X(:,1), X(:,2));

function [C, X] = CreateMesh(xMin, xMax, yMin, yMax, dx, dy)

ix = floor((xMax-xMin)/dx);
jx = floor((yMax-yMin)/dy);
dx = (xMax-xMin)/ix;
dy = (yMax-yMin)/jx;

ii = 1;

for j = 1:jx+1
    for i = 1:ix+1
        X(ii,:) = [ xMin+dx*(i-1), yMin+dy*(j-1)];
        ii = ii+1;
    end
end

C = delaunay(X(:,1), X(:,2));
return;

nL = ix+1;
C = [];
for j = 1:jx
    for i = 1:ix
        t = (j-1)*(ix+1) + i;
        C = [C; t, t+1, t+1+nL; t, t+nL, t+nL];
    end
end




function [Nodes, Elements] = ConvertQuadratic(Nodes, Elements)

nElem = size(Elements,1);
nNodes = size(Nodes,1);

nThis = nNodes;
Elements = [Elements, zeros(nElem, 3)];
for el = 1:nElem
    Cel = Elements(el,1:3);
    Xel = Nodes(Cel,:);
    cNew = [0,0,0];
    for nn = 1:3
        if ( nn == 1)
            this = 1:2;
        elseif ( nn == 2)
            this = 2:3;
        elseif (nn == 3)
            this = [3,1];
        else
            clear this;
        end
        xNew = mean(Xel(this,:));
        index = find( Nodes(:,1) == xNew(1) & Nodes(:,2) == xNew(2));
        if ( isempty(index))
            Nodes = [Nodes; xNew];
            nThis = nThis+1;
            cNew(nn) = nThis;
        else
            cNew(nn) = index;
        end
    end
    Elements(el,4:6) = cNew;
end








