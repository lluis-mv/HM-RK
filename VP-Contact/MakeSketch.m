function MakeSketch(eSize)

if ( nargin == 0)
    eSize = 0.04;
end


model = createpde(1);
R1 = [3,4, 0, 1, 1, 0,  0, 0, -0.25, -0.25]';

g = decsg(R1);
geometryFromEdges(model, g);

model1 = createpde(1);
geometryFromEdges(model1, g);



mesh = generateMesh(model, 'Hmax', eSize);
Nodes = mesh.Nodes';
Elements = mesh.Elements';



mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
Nodes1 = mesh1.Nodes';
Elements1 = mesh1.Elements';

figure(1); clf;
triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')


axis off; axis equal
hold on



B = 2;
x = [0, B/2, B/2, 0];
y = [0,0,B/8,B/8];
pp = polyshape(x, y);
plot(pp, 'facecolor', [0.2, 0.2, 0.2], 'linewidth', 2)
hold on


% plot( [0.01, 0.99], [0.3, 0.3], 'k')
% plot( 0.01, 0.3,  '<k')
% plot( 0.99, 0.3,  '>k')
% tt = text(0.5, 0.35, '$B/2$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');
tt = text(0.5, 0.30, 'Load $q$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');




plot( [-0.1, -0.1], [-0.05 -0.2], 'k')
plot( -0.1, -0.2,  'vk')
plot(-0.1, -0.05,  '^k')
tt = text(-0.05, -0.15, '$B/8$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');

plot( [0.01, 0.99], [-0.3, -0.30], 'k')
plot( 0.01, -0.3,  '<k')
plot( 0.99, -0.3,  '>k')
tt = text(0.5, -0.35, '$B/2$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');


tt = text(-0.2, -0.15, 'Fixed $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', -90)


tt = text(1.1, -0.12, 'Fixed $p$\fontsize{6}{0}\selectfont$_\textnormal{w}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', 90)


tt = text(0.48, -0.44, 'Fixed  $u$\fontsize{6}{0}\selectfont$_\textnormal{v}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');

tt = text(2.7, 0.2, 'Fixed  $p$\fontsize{6}{0}\selectfont$_\textnormal{w}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');


plot(0,0, 'sr', 'linewidth', 1,'MarkerFaceColor', 'r')




print('SketchContact', '-dpdf')
