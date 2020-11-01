

syms alfa real
syms beta real

N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
    alfa*(2*beta-1);
    beta*(2*beta-1);
    4*alfa*(1-alfa-beta);
    4*alfa*beta;
    4*beta*(1-alfa-beta)];

syms a real
syms b real
syms c real
a = sym('aa', [6,1], 'real')


% Na = (a + b*alfa + c*beta)*ones(6,1);

Na = a;

term = Na*(N-Na)'

wA = [2/3,1/6,1/6];
wB = [1/6,1/6,2/3];
w = [1/6,1/6,1/6];

% num integrate
in = 0;
for i = 1:3
    aux = subs(term, alfa, wA(i));
    in = in + subs(aux, beta, wB(i))*w(i);
end