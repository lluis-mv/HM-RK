nSteps = 10.^linspace(1, 5, 1000);

M = 1; 
k = 1;
dt = 0.5./nSteps;

h = 0.05;

AlphaStab = 6*dt*k/h^2;

figure(292)
semilogx(nSteps, AlphaStab, 'r', 'linewidth', 2.5);
hold on

figure(293)
loglog(nSteps, AlphaStab, 'r', 'linewidth', 2.5);
hold on

AlphaStab = 6*dt*k/h^2.*( 1 - exp(-(6*dt*k/h^2).^40 ) );


figure(292)
semilogx(nSteps, AlphaStab, 'b');


figure(293)
loglog(nSteps, AlphaStab, 'b');



figure(292)
hold off
figure(293)
hold off