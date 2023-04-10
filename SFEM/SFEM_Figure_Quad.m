function SFEM_Figure()
addpath('../Sources')
n = 4;

[X, C] = CreateFEM(n);

for n = 1:size(X,1)
    if ( X(n,1) == max(X(:,1)))
        continue
    end
    if ( X(n,1) == min(X(:,1)))
        continue
    end
    if ( X(n,2) == max(X(:,2)))
        continue
    end
    if ( X(n,2) == min(X(:,2)))
        continue
    end
    X(n,:) = X(n,:) + 0.08*rand(1,2);
end

% First part. compute the eigenvalues
figure(1);
clf;



PlotTheMesh(C, X);





[GPNodal, GPElem] = CreateInformation(C, X(:,1)', X(:,2)');


figure(1)

PlotPatches(C, X(:,1)', X(:,2)');
axis off
hold off



figure(2); clf

x = X(:,1)';
y = X(:,2)';

dist = 1000;
theNode = 1;
for i = 1:length(x)
    d = norm( [x(i)-0.5, y(i)-0.5]);
    if ( d<dist)
        dist = d;
        theNode = i;
    end
end

PlotShapeFunction(C, x, y, theNode);
hold on
colorbar
PlotTheMesh(C, X)
axis equal
hold on
axis off;
hold off



figure(3); clf;
PlotGradient(C, x, y, theNode, GPElem);
hold on
colorbar
PlotTheMesh(C, X);
axis equal
hold on
axis off;
hold off




figure(4); clf;
PlotTheMesh(C, X);
hold on
PlotGradientNodal(C, x, y, theNode, GPNodal);
hold on
colorbar
PlotTheMesh(C, X);
axis equal
hold on
axis off;
hold off

for i = 1:4
    figure(i); set(gca, 'FontSize', 15)
end

figure(3);
cc = caxis;
m = max(abs(cc));
cc = [-m, m];
for iii = [3,4]
    figure(iii)
    colormap jet
    caxis(cc);
    colorbar
    drawnow
end



fig = figure(1); drawnow; 
exportgraphics(fig,'QFEM-MeshPatch.pdf', 'BackgroundColor', 'none','ContentType','vector');
fig = figure(2);  drawnow; 
exportgraphics(fig,'QFEM-N.pdf', 'BackgroundColor', 'none','ContentType','vector');
fig = figure(3);  drawnow;
exportgraphics(fig,'QFEM-gradN.pdf', 'BackgroundColor', 'none','ContentType','vector');
fig = figure(4);  drawnow;
exportgraphics(fig,'QFEM-gradS.pdf', 'BackgroundColor', 'none','ContentType','vector');



function PlotShapeFunction(C, x, y, theNode)

nElem = size(C,1)

xP = [];
yP = [];
c = [];
for elem = 1:nElem
    Celem = C(elem,:);

    if ( any(Celem == theNode) )
        index = find(Celem == theNode);

        [al, be] = meshgrid( [-1:0.5:1], [-1:0.5:1]);
        res = nan*al;
        xx = nan*al;
        yy = xx;


        Xel = x(Celem);
        Yel = y(Celem);
        for i = 1:size(al,1)
            for j = 1:size(al,2)
                alfa = al(i,j);
                beta = be(i,j);
                NsmallP =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)];
                xx(i,j) = NsmallP'*Xel';
                yy(i,j) = NsmallP'*Yel';
                res(i,j) = NsmallP(index);
            end
        end
        surf(xx,yy,res,'FaceColor', 'interp' , 'EdgeColor', 'interp')
        view(0,90)
        hold on





    end

end





function PlotGradient(C, x, y, theNode, GPElem)

nElem = size(C,1)

xP = [];
yP = [];
c = [];
for elem = 1:nElem
    Celem = C(elem,:);

    if ( any(Celem == theNode) )
        index = find(Celem == theNode);

        [al, be] = meshgrid( [-1:0.5:1], [-1:0.5:1]);
        res = nan*al;
        xx = nan*al;
        yy = xx;


        Xel = x(Celem);
        Yel = y(Celem);
        X = [Xel', Yel'];
        for i = 1:size(al,1)
            for j = 1:size(al,2)
                alfa = al(i,j);
                beta = be(i,j);
                NsmallP =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)];
                xx(i,j) = NsmallP'*Xel';
                yy(i,j) = NsmallP'*Yel';

                Nsmall_chi = [   beta/4 - 1/4,   alfa/4 - 1/4;
                    1/4 - beta/4, - alfa/4 - 1/4;
                    beta/4 + 1/4,   alfa/4 + 1/4;
                    - beta/4 - 1/4,   1/4 - alfa/4];
                J = Nsmall_chi'*X;
                dN_dX = inv(J)*Nsmall_chi';

                res(i,j) = dN_dX(1, index);
            end
        end
        surf(xx,yy,res,'FaceColor', 'interp' , 'EdgeColor', 'interp')
        view(0,90)
        hold on





    end

end







function PlotGradientNodal(C, x, y, theNode, GPNodal)

nNodes = length(x);

xP = [];
yP = [];
c = [];

for node = 1:nNodes
    if (any(GPNodal(node).NeigNodes == theNode))
        index = find( GPNodal(node).NeigNodes == theNode);
        value = GPNodal(node).dN_dX(1,index);

        for elem = [GPNodal(node).NeigElement]'

            Celem = C(elem,:);
            index = find(Celem == node);
            if ( index == 1)
                i1 = Celem(4); i2 = Celem(2);
            elseif ( index == 2)
                i1 = Celem(1); i2 = Celem(3);
            elseif ( index == 3)
                i1 = Celem(2); i2 = Celem(4);
            else
                i1 = Celem(3); i2 = Celem(1);
            end

            x1 = [x(node), mean( [x(i2), x(node)]), mean(x(Celem)), mean([x(i1), x(node)])];
            y1 = [y(node), mean( [y(i2), y(node)]), mean(y(Celem)), mean([y(i1), y(node)])];
            xP = [xP, x1'];
            yP = [yP, y1'];
            c = [c, value];


        end

        patch(xP, yP, c)
    end
end



function [] = PlotPatches(C, x, y)

nElem = size(C,1);

for elem = 1:nElem


    Celem = C(elem,:);
    xC = mean(x(Celem));
    yC = mean(y(Celem));

    x1 = mean(x(Celem(1:2)));
    x2 = mean(x(Celem([2,3])));
    x3 = mean(x(Celem(3:4)));
    x4 = mean(x(Celem([4,1])));

    y1 = mean(y(Celem(1:2)));
    y2 = mean(y(Celem(2:3)));
    y3 = mean(y(Celem(3:4)));
    y4 = mean(y(Celem([4,1])));


    plot( [x1, xC], [y1, yC], 'b-.')
    plot( [x2, xC], [y2, yC], 'b-.')
    plot( [x3, xC], [y3, yC], 'b-.')
    plot( [x4, xC], [y4, yC], 'b-.')




end


function [GPNodal, GPElem] = CreateInformation(C, x, y)

CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-2;
CP.Elastic = false;
CP.MCC = true;

Elements = C;
ElementType = 'Q4Q4';
Nodes = [x;y]';

[GPElem] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPElem] = InitializeConstitutiveLaw(CP, GPElem);
GPElem = CalculateSmoothingPathAreas( Nodes, Elements, GPElem)
[GPNodal] = ConstructNodalIntegrationPointsQuad(CP, Nodes, Elements, GPElem);






function [Nodes, Elements] = CreateFEM( nx)

ny = floor( nx);

dx = 1/nx;
dy = 1/ny;

nx = nx+1;
ny = ny+1;

Elements = [];
for j = 1:ny-1
    for i = 1:nx-1

        Elements = [ Elements;
            (j-1)*nx+i, (j-1)*nx+i+1, (j)*nx+i+1, (j)*nx+i];
    end
end

Nodes = [];
for j = 1:ny
    for i = 1:nx
        Nodes = [ Nodes;
            (i-1)*dx, (j-1)*dy];
    end
end

function PlotTheMesh(C, X)
ii = [1,2,3,4,1];
for elem = 1:size(C,1)
    Celem = C(elem,:);
    Xe = X(Celem,:);
    plot3(Xe(ii,1), Xe(ii,2), 50*Xe(ii,1), 'k', 'linewidth', 3)
    hold on
end
view(0,90)
drawnow
axis equal
axis off
