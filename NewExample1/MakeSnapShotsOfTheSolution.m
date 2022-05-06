function [] = MakeSnapShotsOfTheSolution()

addpath('../Sources')
% 1. Define the problem




CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




for i = 1:2
    figure(i);
    clf;
     newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    colororder(newcolors)
end


for dt = [1E-4, 1E-2, 1E-1, 1, 10]
[Nodes, Elements] = ReadTheMesh('ThisMesh.msh');
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'Q8Q4');
U = zeros(3*size(Nodes,1), 1);
[Xa] = ComputeAnalyticalSolution(Nodes, Elements,'Q8Q4', dt, CP, GPInfo, U);





figure(1)
index = find(Nodes(:,1) == 0);
Uv = Xa(3*(index-1)+2);
y = Nodes(index,2);
[ya, index] = sort(y);
Uv = Uv(index);
% [Uv, ya] = CorrectInterpolation(Uv, ya);
plot(Uv, ya, '', 'linewidth', 2, 'DisplayName', 'Reference')
hold on
drawnow


figure(2)
index = find(Nodes(:,1) == 0);
WP = Xa(3*(index-1)+3);
y = Nodes(index,2);
[ya, index] = sort(y);
WPa = WP(index);
% [WPa, ya] = CorrectInterpolation(WPa, ya);
plot(WPa, ya, '', 'linewidth', 2, 'DisplayName', 'Reference')
hold on
drawnow
end

figure(1)
xlabel('$u_v$ (m)', 'interpreter', 'latex')
ylabel('$y$ (m)', 'interpreter', 'latex')
set(gca, 'FontSize', 15)
figure(2)
xlabel('$p_w$ (kPa)', 'interpreter', 'latex')
ylabel('$y$ (m)', 'interpreter', 'latex')
set(gca, 'FontSize', 15)

for fig = 1:2
    figure(fig)
    legend('$t = 10^{-4}$ s', '$t = 10^{-2}$ s', '$t = 10^{-1}$ s', '$t = 1$ s', '$t = 10$ s', 'interpreter', 'latex', 'location', 'best')
    drawnow;
    print(['ReferenceSolution-', num2str(fig)], '-dpdf')
end






function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 1, t, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, Ca, Ka);

Aa = Ca\(Ka);

% [vectors, values] = eig(full(Aa), 'nobalance');
[vectors, values] = eig(full(Aa), eye(size(Aa)), 'qz');
Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end

Xa = real(Xa);


