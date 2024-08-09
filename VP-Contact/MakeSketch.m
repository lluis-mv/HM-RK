function MakeSketch(XFILE)
addpath('../Sources')
if ( nargin == 0)
    XFILE = 'Mesh.msh';
end


[Nodes, Elements] = ReadTheMesh('Mesh.msh');
figure(1); clf;
triplot(Elements(:,1:3), Nodes(:,1), Nodes(:,2), 'k')


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
