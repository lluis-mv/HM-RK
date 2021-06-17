

clf;
poisson = [-1:0.01:0.5];
M = 1;

tauLMV = 0.0*poisson; tauLI = 0.0*poisson; tauBorja = 0.0*poisson;
for i = 1:length(poisson)
    nu = poisson(i);
    tauLMV(i) = 2/M;
    tauLI(i) = 3/M;
    G = M * ( 1- 2 * nu) / 2 / (1-nu);
    tauBorja(i) = 1/(2*G);
end

plot(poisson, tauLMV)
hold on
plot(poisson, tauLI)
plot(poisson, tauBorja)
legend('LMV', 'LI', 'BORJA')