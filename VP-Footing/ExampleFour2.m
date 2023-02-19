function [] = ExampleFour2()

load('UndrainedData.mat')

% eSizeAxis = ESIZE(1:i);
figure(99); clf
semilogx(eSizeAxis, Qnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(eSizeAxis, Qlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(eSizeAxis, Qmixed, 'cs-.', 'DisplayName', 'T3T3T3')
plot(eSizeAxis, Qquad, 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('$h_e$ (m)', 'interpreter', 'latex')
ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
ylim([32,37])


figure(104); clf
loglog(eSizeAxis, abs(Qnodal-Qquad(end)), 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(eSizeAxis, abs(Qlinear-Qquad(end)), 'gs-.', 'DisplayName', 'T3T3')
plot(eSizeAxis, abs(Qmixed-Qquad(end)), 'cs-.', 'DisplayName', 'T3T3T3')
plot(eSizeAxis, abs(Qquad-Qquad(end)), 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('$h_e$ (m)', 'interpreter', 'latex')
ylabel('Error (kPa)', 'interpreter', 'latex')


figure(100); clf
semilogx(eSizeAxis, PWnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(eSizeAxis, PWlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(eSizeAxis, PWmixed, 'cs-.', 'DisplayName', 'T3T3T3')
plot(eSizeAxis, PWquad, 'bs-.', 'DisplayName', 'T6T3')

xlabel('$h_e$ (m)', 'interpreter', 'latex')
ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
drawnow
ylim([18,23])



figure(101); clf
loglog(nDofs, TIMEnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, TIMElinear, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsmixed, TIMEmixed, 'c*-.', 'DisplayName', 'T3T3T3')
plot(nDofsquad, TIMEquad, 'b*-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('Computational cost (s)', 'interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


figure(102); clf
semilogx(nDofs, TIMEnodal./nDofs, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, TIMElinear./nDofs, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsmixed, TIMEmixed./nDofsmixed, 'c*-.', 'DisplayName', 'T3T3T3')
plot(nDofsquad, TIMEquad./nDofsquad, 'b*-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('Computational cost / $n_{dofs}$ (s)', 'interpreter', 'latex')
set(gca, 'XScale', 'log')
drawnow


figure(103); clf
semilogx(nDofs, nZeronodal./nDofs, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, nZerolinear./nDofs, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsmixed, nZeromixed./nDofsmixed, 'c*-.', 'DisplayName', 'T3T3T3')
plot(nDofsquad, nZeroquad./nDofsquad, 'b*-.', 'DisplayName', 'T6T3')

xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('$n_{nz}$/ $n_{dofs}$ (s)', 'interpreter', 'latex')
drawnow


figure(200); clf
loglog(TIMEnodal, abs(Qnodal-Qquad(end)), 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(TIMElinear, abs(Qlinear-Qquad(end)), 'gs-.', 'DisplayName', 'T3T3')
plot(TIMEmixed, abs(Qmixed-Qquad(end)), 'cs-.', 'DisplayName', 'T3T3T3')
plot(TIMEquad, abs(Qquad-Qquad(end)), 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Computational cost (s)', 'interpreter', 'latex')
ylabel('Error (kPa)', 'interpreter', 'latex')
drawnow
% ylim([32,37])

figure(201); clf
loglog(TIMEnodal, abs(Qnodal-0*Qquad(end)), 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(TIMElinear, abs(Qlinear-0*Qquad(end)), 'gs-.', 'DisplayName', 'T3T3')
plot(TIMEmixed, abs(Qmixed-0*Qquad(end)), 'cs-.', 'DisplayName', 'T3T3T3')
plot(TIMEquad, abs(Qquad-0*Qquad(end)), 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Computational cost (s)', 'interpreter', 'latex')
ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
drawnow
ylim([32,37])


figure(202); clf
loglog(TIMEnodal, abs(PWnodal-PWquad(end)), 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(TIMElinear, abs(PWlinear-PWquad(end)), 'gs-.', 'DisplayName', 'T3T3')
plot(TIMEmixed, abs(PWmixed-PWquad(end)), 'cs-.', 'DisplayName', 'T3T3T3')
plot(TIMEquad, abs(PWquad-PWquad(end)), 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Computational cost (s)', 'interpreter', 'latex')
ylabel('Error (kPa)', 'interpreter', 'latex')
drawnow
% ylim([32,37])

figure(203); clf
loglog(TIMEnodal, abs(PWnodal-0*PWquad(end)), 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(TIMElinear, abs(PWlinear-0*PWquad(end)), 'gs-.', 'DisplayName', 'T3T3')
plot(TIMEmixed, abs(PWmixed-0*PWquad(end)), 'cs-.', 'DisplayName', 'T3T3T3')
plot(TIMEquad, abs(PWquad-0*PWquad(end)), 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Computational cost (s)', 'interpreter', 'latex')
ylabel('Water pressure (kPa)', 'interpreter', 'latex')
drawnow
ylim([18,23])



for jj = [99, 100, 101, 102, 103, 104, 200:203]
    figure(jj)
    set(gca, 'FontSize', 15)
    ll = legend();
    set(ll, 'location', 'best', 'interpreter', 'latex')
end

figure(99); drawnow; pause(1); print('Footing-Load-1', '-dpdf')
figure(104); drawnow; pause(1); print('Footing-Convergence-1', '-dpdf')
figure(100); drawnow; pause(1); print('Footing-WP-1', '-dpdf')
figure(101); drawnow; pause(1); print('Footing-Cost-1', '-dpdf')
figure(102); drawnow; pause(1); print('Footing-Velocity-1', '-dpdf')
figure(103); drawnow; pause(1); print('Footing-Nnz-1', '-dpdf')
figure(200); drawnow; pause(1); print('Footing-Efficiency-1', '-dpdf')
figure(201); drawnow; pause(1); print('Footing-Efficiency-2', '-dpdf')
figure(202); drawnow; pause(1); print('Footing-Efficiency-3', '-dpdf')
figure(203); drawnow; pause(1); print('Footing-Efficiency-4', '-dpdf')


close all; clear all;

load('TimeUndrainedData.mat')

figure(99); clf
plot(NSTEPS, Qnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(NSTEPS, Qlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(NSTEPS, Qmixed, 'cs-.', 'DisplayName', 'T3T3T3')
plot(NSTEPS, Qquad, 'bs-.', 'DisplayName', 'T6T3')
drawnow
xlabel('$n_{steps}$', 'interpreter', 'latex')
ylabel('Footing resistance (kPa)', 'interpreter', 'latex')
ylim([32,37])
set(gca, 'XScale', 'log')


figure(100); clf
plot(NSTEPS, PWnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(NSTEPS, PWlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(NSTEPS, PWmixed, 'cs-.', 'DisplayName', 'T3T3T3')
plot(NSTEPS, PWquad, 'bs-.', 'DisplayName', 'T6T3')
xlabel('$n_{steps}$', 'interpreter', 'latex')
ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
ylim([18,23])
drawnow
set(gca, 'XScale', 'log')

for jj = [99, 100, ]
    figure(jj)
    set(gca, 'FontSize', 15)
    ll = legend();
    set(ll, 'location', 'best', 'interpreter', 'latex')
end

figure(99); drawnow; pause(1); print('Footing-Load-2', '-dpdf')
figure(100); drawnow; pause(1); print('Footing-WP-2', '-dpdf')
