
syms alfa real
syms beta real
syms h real


A = int(term, beta, [0, 1-alfa]);
A = int(A, alfa, [0, 1])

N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
    alfa*(2*alfa-1);
    beta*(2*beta-1);
    4*alfa*(1-alfa-beta);
    4*alfa*beta;
    4*beta*(1-alfa-beta)];



syms a real;
syms b real
syms c real
Na = (a+b*alfa+c*beta)*ones(6,1);


term = Na*(N'-Na')

A = int(term, beta, [0, 1-alfa]);
A = int(A, alfa, [0, 1])
%
% a2 = solve(A(1,1), a);
% a2 = a2(1);
%
% A2 = subs(A, a, a2)
% b2 = solve(A(2,1), b)
% b2 = b2(1),
%
% A3 = subs(A2, b, b2);
% A3 = simplify(A3)