function SFEM_Figure()
addpath('../')
n = 5;

x = meshgrid(linspace(0,1,n));
[x,y] = meshgrid(linspace(0,1,n));
x = x(1:size(x,1)*size(x,2));
y = y(1:size(x,1)*size(x,2));
x = x+0.02*rand(1,length(x));
y = y+0.02*rand(1,length(x));



C = delaunay(x, y);


[GPNodal, GPElem] = CreateInformation(C, x, y);


figure(1)
triplot(C, x,y, 'k', 'linewidth', 3)
axis equal
hold on
PlotPatches(C, x, y)
axis off
hold off



figure(2); clf



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
triplot(C, x,y, 'k', 'linewidth', 3)
axis equal
hold on
axis off;
hold off



figure(3); clf;
PlotGradient(C, x, y, theNode, GPElem);
hold on
colorbar
triplot(C, x,y, 'k', 'linewidth', 3)
axis equal
hold on
axis off;
hold off




figure(4); clf;
triplot(C, x,y, 'k', 'linewidth', 3)
hold on
PlotGradientNodal(C, x, y, theNode, GPNodal);
hold on
colorbar
triplot(C, x,y, 'k', 'linewidth', 3)
axis equal
hold on
axis off;
hold off

function PlotShapeFunction(C, x, y, theNode)

nElem = size(C,1)

xP = [];
yP = [];
c = [];
for elem = 1:nElem
    Celem = C(elem,:);
    
    if ( any(Celem == theNode) )
        xP = [xP, x(Celem)'];
        yP = [yP, y(Celem)'];
        
        index = find(Celem == theNode);
        v = zeros(3,1); v(index) = 1;
        c = [c, v];
        
        
    end

end

patch(xP, yP, c, 'FaceColor', 'interp')



function PlotGradient(C, x, y, theNode, GPElem)

nElem = size(C,1)

xP = [];
yP = [];
c = [];
for elem = 1:nElem
    Celem = C(elem,:);
    
    if ( any(Celem == theNode) )
        xP = [xP, x(Celem)'];
        yP = [yP, y(Celem)'];
        
        
        index = find(Celem == theNode);
        v = GPElem(elem).dN_dX(1,index);
        c = [c, v];
        
        
    end

end

patch(xP, yP, c)


function PlotGradientNodal(C, x, y, theNode, GPNodal)

nNodes = length(x);

xP = [];
yP = [];
c = [];

for node = 1:nNodes
    if (any(GPNodal(node).NeigNodes == theNode))
        index = find( GPNodal(node).NeigNodes == theNode);
        value = GPNodal(node).dN_dX(1,index);
% %         
% %         % now I have to create something to plot the patch value....
% %         u = []; v = [];
% %         for nNode = [GPNodal(node).NeigNodes]'
% %             u = [u, mean([x(node),x(nNode)])]
% %             v = [v, mean([y(node),y(nNode)])]
% %         end
% %         
% %         for elem = [GPNodal(node).NeigElement]'
% %             Celem = C(elem,:);
% %             u = [u, mean( x(Celem))];
% %             v = [v, mean( y(Celem))];
% %         end
% %         
% %         tri = delaunay(u, v);
% %         
% %         for m = 1:size(tri,1)
% %             xP = [xP, u(tri(m,:))']
% %             yP = [yP, v(tri(m,:))']
% %             c = [c, value];
% %         end
        
        for elem = [GPNodal(node).NeigElement]'
            
            Celem = C(elem,:)
            index = find(Celem == node);
            if ( index == 1)
                i1 = Celem(2); i2 = Celem(3);
            elseif ( index == 2)
                i1 = Celem(1); i2 = Celem(3);
            else
                i1 = Celem(1); i2 = Celem(2);
            end
            
            x1 = [x(node), mean( [x(i1), x(node)]), mean(x(Celem)), mean([x(i2), x(node)])];
            y1 = [y(node), mean( [y(i1), y(node)]), mean(y(Celem)), mean([y(i2), y(node)])];
            xP = [xP, x1'];
            yP = [yP, y1'];
            c = [c, value];
            
            
        end

        patch(xP, yP, c)        
    end
end



function [] = PlotPatches(C, x, y)

nElem = size(C,1)

for elem = 1:nElem
    Celem = C(elem,:);
    xC = mean(x(Celem));
    yC = mean(y(Celem));
    
    x1 = mean(x(Celem(1:2)));
    x2 = mean(x(Celem([1,3])));
    x3 = mean(x(Celem(2:3)));
    
    y1 = mean(y(Celem(1:2)));
    y2 = mean(y(Celem([1,3])));
    y3 = mean(y(Celem(2:3)));
    
    plot( [x1, xC], [y1, yC], 'b-.')
    plot( [x2, xC], [y2, yC], 'b-.')
    plot( [x3, xC], [y3, yC], 'b-.')
    
    
    
    
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
ElementType = 'T3T3';
Nodes = [x;y]';

[GPElem] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPElem] = InitializeConstitutiveLaw(CP, GPElem);
[GPNodal] = ConstructNodalIntegrationPoints(CP, Nodes, Elements, GPElem);


