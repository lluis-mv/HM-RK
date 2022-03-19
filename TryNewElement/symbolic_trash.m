syms alfa real
syms beta real

Nsmall =  [ -1/4*(1-alfa)*(1-beta)*(1+alfa+beta);
    -1/4*(1+alfa)*(1-beta)*(1-alfa+beta);
    -1/4*(1+alfa)*(1+beta)*(1-alfa-beta);
    -1/4*(1-alfa)*(1+beta)*(1+alfa-beta);
    1/2*(1-alfa^2)*(1-beta);
    1/2*(1+alfa)*(1-beta^2);
    1/2*(1-alfa^2)*(1+beta);
    1/2*(1-alfa)*(1-beta^2);
    ];

Nsmall_chi = simplify([diff(Nsmall, alfa), diff(Nsmall, beta)])

X = [0,0;
    1, 0;
    1, 1;
    0, 1;
    0.5, 0;
    1, 0.5;
    0.5, 1;
    0, 0.5];

J = simplify(Nsmall_chi'*X)

deris = inv(J)*Nsmall_chi'