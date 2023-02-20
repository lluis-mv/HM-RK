function [] = CompareElements()

new_figure(212); hold off; clf;
new_figure(214); hold off; clf;
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = false;
CP.MCC = true;

% Modified Cam clay parameters
CP.kappa = 0.01;
CP.lambda = 0.1;
CP.M_MCC = 1;
CP.nu = 0.3;
CP.u = 20;



CP.RK = 1;


LinearElastic = false;
if ( LinearElastic)
    CP.Elastic = true;
    CP.MCC = false;
end



[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);


figure(88); clf
PlotMesh(NodesQ, ElementsQ, false)
drawnow
axis equal
axis off
print('MeshFooting-Q', '-dpdf')
figure(89); clf
PlotMesh(NodesT, ElementsT, false)
drawnow
axis equal
axis off
print('MeshFooting-T', '-dpdf')
figure(100); clf; 



B = 1;

x = [0, B/2, B/2, 0];
y = [0,0,B/4,B/4];
pp = polyshape(x, y);
plot(pp, 'facecolor', [0.2, 0.2, 0.2], 'linewidth', 2)
hold on

% plot( [-0.1, -0.1], [0.05 0.45], 'k')
% plot( -0.1, 0.05,  'vk')
% plot(-0.1, 0.45,  '^k')
%tt = text(-0.3, 0.25, '$B/2$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');

plot( [0.05, 0.45], [0.35, 0.35], 'k')
plot( 0.05, 0.35,  '<k')
plot( 0.45, 0.35,  '>k')
tt = text(0.25, 0.45, '$B/2$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');


B = 0.4
x = 5*B*[0, 1, 1, 0];
y = 3.8*B*[0,0,-1, -1];
pp = polyshape(x, y);
plot(pp, 'facecolor', 'w', 'linewidth', 2)


plot( [-0.1, -0.1], [-0.05 -1.5], 'k')
plot( -0.1, -1.5,  'vk')
plot(-0.1, -0.05,  '^k')
tt = text(-0.3, -0.75, '$2 B$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');

plot( [0.05, 1.95], [-1.6, -1.6], 'k')
plot( 0.05, -1.6,  '<k')
plot( 1.95, -1.6,  '>k')
tt = text(1, -1.7, '$2\,B$', 'interpreter', 'latex','fontsize', 12, 'HorizontalAlignment', 'Center');


tt = text(0.1, -0.75, 'Fixed $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', -90)


tt = text(1.9, -0.75, 'Fixed $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');
set(tt, 'rotation', 90)


tt = text(1, -1.4, 'Fixed  $u$\fontsize{6}{0}\selectfont$_\textnormal{h}$\fontsize{10}{0}\selectfont, $u$\fontsize{6}{0}\selectfont$_\textnormal{v}$ \fontsize{10}{0}\selectfont and $p$\fontsize{6}{0}\selectfont$_\textnormal{w}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');

tt = text(1.2, -0.2, 'Fixed  $p$\fontsize{6}{0}\selectfont$_\textnormal{w}$', 'interpreter', 'latex','fontsize', 14, 'HorizontalAlignment', 'Center');


plot(0,0, 'sr', 'linewidth', 1,'MarkerFaceColor', 'r')

axis equal
axis off


print('Sketch1', '-dpdf')




RK = 1;


nSteps = 100;
dt = 1/nSteps;


%%%%%%%%%%%%%%%%%%%%%%%%% T6T3 

[U1, GPInfo1, rrr1,  information1] = ComputeNLProblem(NodesT, ElementsT, CP, dt, nSteps, 'T6T3', RK, 1, false);



new_figure(1)
PlotNodal(NodesT, ElementsT, U1(3:3:end))
drawnow;
colorbar; 




FF = [information1.F];
new_figure(212)
plot( FF(3:3:end), FF(1:3:end), 'r', 'linewidth', 2 , 'DisplayName', ['T6T3'])
hold on
drawnow

new_figure(214)
plot( FF(3:3:end), FF(2:3:end), 'r', 'linewidth', 2, 'DisplayName', ['T6T3'])
hold on
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%% Q8Q4





[U2, GPInfo, rrr,  information2] = ComputeNLProblem(NodesQ, ElementsQ, CP, dt, nSteps, 'Q8Q4', RK, 1, false);


new_figure(2)
PlotNodal(NodesQ, ElementsQ, U2(3:3:end))
drawnow;
colorbar; 

FF = [information2.F];
new_figure(212)
plot( FF(3:3:end), FF(1:3:end), 'g', 'linewidth', 2,'DisplayName', ['Q8Q4'])
hold on
drawnow

new_figure(214)
plot( FF(3:3:end), FF(2:3:end), 'g', 'linewidth', 2,'DisplayName', ['Q8Q4'])
hold on
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%% T6T3  Implicit


[Uimp, GPInfo, rrr,  informationI] = ComputeImplicitNonLinearProblem(NodesT, ElementsT, CP, dt, nSteps, 'T6T3');

new_figure(3)
PlotNodal(NodesT, ElementsT, Uimp(3:3:end))
drawnow;
colorbar; 


FF = [informationI.F];
new_figure(212)
plot( FF(3:3:end), FF(1:3:end), 'k:', 'linewidth', 2, 'DisplayName',  ['T6T3 Implicit'])
hold on
drawnow

new_figure(214)
plot( FF(3:3:end), FF(2:3:end), 'k:', 'linewidth', 2, 'DisplayName', ['T6T3 Implicit'])
hold on
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%% Q8Q4  Implicit


[Uimp2, GPInfo, rrr,  informationI] = ComputeImplicitNonLinearProblem(NodesQ, ElementsQ, CP, dt, nSteps, 'Q8Q4');

new_figure(4)
PlotNodal(NodesQ, ElementsQ, Uimp2(3:3:end) )
drawnow;
colorbar; 


FF = [informationI.F];
new_figure(212)
plot( FF(3:3:end), FF(1:3:end), 'b:', 'linewidth', 2, 'DisplayName',  ['Q8Q4 Implicit'])
hold on
drawnow

new_figure(214)
plot( FF(3:3:end), FF(2:3:end), 'b:', 'linewidth', 2, 'DisplayName', ['Q8Q4 Implicit'])
hold on
drawnow






% writting everything down
    
maxx = 0;
for i = [1:4]
    new_figure(i)
    xx = caxis();
    maxx = maxx+xx(2);
end

for i = [1:4]
    new_figure(i)
    drawnow
    caxis([0, maxx/4])
    axis equal
    xlim([0,4])
    ylim([-4,0])
    axis off
    pbaspect([1 1 10])
    drawnow
    MyPrint(['WaterPressurePlastic-', num2str(i)], '-dpdf')
end

new_figure(212)
drawnow

xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Footing pressure (kN)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
xlim([0, 0.05])
drawnow
MyPrint('Plastic-RR', '-dpdf')

new_figure(214)
drawnow
xlabel('Footing indentation (m)', 'interpreter', 'latex')
ylabel('Water pressure (kPa)', 'interpreter', 'latex')
set(gca, 'FontSize', 13)
legend('location', 'best', 'interpreter', 'latex')
xlim([0, 0.05])
drawnow
MyPrint('Plastic-WP', '-dpdf')


function [] = MyPrint( NAME, FORMAT)

drawnow;
print(NAME, FORMAT)
function [] = new_figure(fig)
fig = figure(fig);
fig.Renderer='Painters';
