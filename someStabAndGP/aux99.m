syms alfa real
syms beta real
syms a real
syms b real


N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
    alfa*(2*alfa-1);
    beta*(2*beta-1);
    4*alfa*(1-alfa-beta);
    4*alfa*beta;
    4*beta*(1-alfa-beta)];


Na = [a,a,a,b,b,b]';

term = Na*(N-Na)'

I = int(term, beta, [0, 1-alfa]);
I = int(I, alfa, [0, 1])




% Numerical integral

wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;



auxK = [Na2, Na2, wa;
Na1, Na2, wa;
Na2, Na1, wa;
Nb2, Nb2, wb ;
Nb1, Nb2, wb ;
Nb2, Nb1, wb ];

wA = auxK(:,1)';
wB = auxK(:,2)';
wW = auxK(:,3)';


In = zeros(6,6);
for i = 1:length(wW)
    
    In = In + wW(i)*subs(subs( term, alfa, wA(i)),beta, wB(i));
end
    
