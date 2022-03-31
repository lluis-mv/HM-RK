close all;

fig1 = figure(1);
fig1.Renderer = 'painters';

[x, y] = meshgrid( 0:0.1:1, 0:0.1:1);


z = rand(size(x));

surf(x,y,z,'FaceColor','interp', 'edgecolor', 'interp')
axis equal
axis off
view(0,90)
print('WithWhiteLines', '-dpdf')


hold on
surf(x,y,z,'FaceColor','interp', 'edgecolor', 'interp')
print('WithoutWhiteLines', '-dpdf')



fig2 = figure(2);
surf(x,y,z,'FaceColor','interp', 'edgecolor', 'interp')
axis equal
axis off
view(0,90)
print('WithoutPainters', '-dpdf')