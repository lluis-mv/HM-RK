

load('UndrainedData.mat')

% eSizeAxis = ESIZE(1:i);
figure(99); clf
semilogx(eSizeAxis, Qnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(eSizeAxis, Qlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(eSizeAxis, Qquad, 'bs-.', 'DisplayName', 'T6T3')
plot(eSizeAxis(NStabilizedNodal), Qnodal(NStabilizedNodal), 'rs','MarkerFaceColor', 'r', 'HandleVisibility', 'off')
plot(eSizeAxis(NStabilizedT3T3), Qlinear(NStabilizedT3T3), 'gs','MarkerFaceColor', 'g', 'HandleVisibility', 'off')
drawnow
xlabel('$h_e$ (m)', 'interpreter', 'latex')
ylabel('Footing resistance (kPa)', 'interpreter', 'latex')

figure(100); clf
semilogx(eSizeAxis, PWnodal, 'rs-.', 'DisplayName', 'NS-T3T3')
hold on
plot(eSizeAxis, PWlinear, 'gs-.', 'DisplayName', 'T3T3')
plot(eSizeAxis, PWquad, 'bs-.', 'DisplayName', 'T6T3')

plot(eSizeAxis(NStabilizedNodal), PWnodal(NStabilizedNodal), 'rs','MarkerFaceColor', 'r', 'HandleVisibility', 'off')
plot(eSizeAxis(NStabilizedT3T3), PWlinear(NStabilizedT3T3), 'gs','MarkerFaceColor', 'g', 'HandleVisibility', 'off')
xlabel('$h_e$ (m)', 'interpreter', 'latex')
ylabel('$p_w$ (kPa)', 'interpreter', 'latex')
drawnow

figure(101); clf
loglog(nDofs, TIMEnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, TIMElinear, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsquad, TIMEquad, 'b*-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('Computational cost (s)', 'interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


figure(102); clf
loglog(nDofs, TIMEnodal./nDofs, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, TIMElinear./nDofs, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsquad, TIMEquad./nDofsquad, 'b*-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('Computational cost / $n_{dofs}$ (s)', 'interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


figure(103); clf
semilogx(nDofs, nZeronodal./nDofs, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(nDofs, nZerolinear./nDofs, 'g*-.', 'DisplayName', 'T3T3')
plot(nDofsquad, nZeroquad./nDofsquad, 'b*-.', 'DisplayName', 'T6T3')

xlabel('Number of dofs, $n_{dofs}$', 'interpreter', 'latex')
ylabel('$n_{nz}$/ $n_{dofs}$ (s)', 'interpreter', 'latex')


figure(105); clf
plot(abs(Qnodal-Qquad(end)), TIMEnodal, 'r*-.', 'DisplayName', 'NS-T3T3')
hold on
plot(abs(Qlinear-Qquad(end)), TIMElinear, 'g*-.', 'DisplayName', 'T3T3')
plot(abs(Qquad-Qquad(end)), TIMEquad, 'b*-.', 'DisplayName', 'T6T3')
drawnow
xlabel('Number of dofs', 'interpreter', 'latex')
ylabel('Computational cost (s)', 'interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


for jj = [99, 100, 101, 102, 103]
    figure(jj)
    set(gca, 'FontSize', 15)
    ll = legend();
    set(ll, 'location', 'best', 'interpreter', 'latex')
end

figure(99); drawnow; pause(1); print('Footing-Load-1', '-dpdf')
figure(100); drawnow; pause(1); print('Footing-WP-1', '-dpdf')
figure(101); drawnow; pause(1); print('Footing-Cost-1', '-dpdf')
figure(102); drawnow; pause(1); print('Footing-Velocity-1', '-dpdf')
