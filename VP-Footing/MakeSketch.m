




model = createpde(1);
R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';

g = decsg(R1);
geometryFromEdges(model, g);

model1 = createpde(1);
geometryFromEdges(model1, g);

eSize = 1/4;

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
y = [0,0,B/4,B/4];
pp = polyshape(x, y);
plot(pp, 'facecolor', [0.2, 0.2, 0.2], 'linewidth', 2)
hold on


plot( [0.05, 0.9], [0.7, 0.7], 'k')
plot( 0.05, 0.7,  '<k')
plot( 0.9, 0.7,  '>k')
tt = text(0.5, 0.9, '$B/2$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');





plot( [-0.1, -0.1], [-0.05 -3.9], 'k')
plot( -0.1, -3.9,  'vk')
plot(-0.1, -0.05,  '^k')
tt = text(-0.3, -2.0, '$2 B$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');

plot( [0.05, 3.95], [-4.2, -4.2], 'k')
plot( 0.05, -4.2,  '<k')
plot( 3.95, -4.2,  '>k')
tt = text(0.5, -4.5, '$2\,B$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');


tt = text(-0.5, -0.75, 'Fixed $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', -90)


tt = text(4.3, -1.75, 'Fixed $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', 90)


tt = text(2.5, -4.5, 'Fixed  $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$\fontsize{10}{0}\selectfont  $\,$ and $u$\fontsize{6}{0}\selectfont$_\textnormal{v}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');

tt = text(2.7, 0.2, 'Fixed  $p$\fontsize{6}{0}\selectfont$_\textnormal{w}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');


plot(0,0, 'sr', 'linewidth', 1,'MarkerFaceColor', 'r')




print('SketchFooting', '-dpdf')
