function auxT6()

a1 = [1/3];
b1 = [1/3];
w1 = [1/2];


a3 = [1/6, 2/3, 1/6];
b3 = [1/6, 1/6, 2/3];
w3 = [1/6, 1/6, 1/6];

a4 = [1/3,3/5,1/5,1/5];
b4 = [1/3,1/5,3/5,1/5];
w4 = [-27/96, 25/96*ones(1,3)];



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

aK = auxK(:,1)';
bK = auxK(:,2)';
wK = auxK(:,3)';


[M1] = Integrate2(a1,b1,w1)

[M3] = Integrate2(a3,b3,w3)

[M4] = Integrate2(a4,b4,w4)
[Mk] = Integrate2(aK,bK,wK)






Ma = IntegrateA();
Ma = eval(Ma)

% M = IntegrateHardMethod()

hola = 1;

function [M] = Integrate2(a,b,w)

M = zeros(6,6);

x0 = [ 2/5, -1/5, -1/5, 3/5, -1/5, 3/5 -3/5, 3/5, 0, 0, 4/5, -4/5 -3/5, 0, 3/5, -4/5, 4/5, 0]';
ap = x0(1:6);
bp = x0(7:12);
cp = x0(13:18);

for i = 1:length(w)
    alfa = a(i);
    beta = b(i);
    N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];
    Na = ap+bp*alfa+cp*beta;
    
    term = N*N'-Na*Na';
    
    M = M + term*w(i);
end





function [M] = IntegrateHardMethod

    x0 = [ 2/5, -1/5, -1/5, 3/5, -1/5, 3/5, ...
    -3/5, 3/5, 0, 0, 4/5, -4/5 ...
    -3/5, 0, 3/5, -4/5, 4/5, 0]';
ap = x0(1:6);
bp = x0(7:12);
cp = x0(13:18);
ns = 10000
i = 0;
M = zeros(6,6);
while (i<ns)
    alfa = rand();
    beta = rand();
    if ( beta > 1-alfa)
        continue;
    end
    N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];
    

Na = ap+bp*alfa+cp*beta;

   term = N*N'-Na*Na';
    M = M + term/ns*0.5;
    i = i+1;
end



function [M] = IntegrateA()

syms alfa real
syms betta real



N =  [ (1 - alfa - betta)*(1-2*alfa-2*betta);
    alfa*(2*alfa-1);
    betta*(2*betta-1);
    4*alfa*(1-alfa-betta);
    4*alfa*betta;
    4*betta*(1-alfa-betta)];

x0 = [ 2/5, -1/5, -1/5, 3/5, -1/5, 3/5, ...
    -3/5, 3/5, 0, 0, 4/5, -4/5 ...
    -3/5, 0, 3/5, -4/5, 4/5, 0]';
ap = x0(1:6);
bp = x0(7:12);
cp = x0(13:18);
Na = ap+bp*alfa+cp*betta;

term =N*(N)'-Na*Na';
M = int(term, betta, [0, 1-alfa]);
M = int(M, alfa, [0, 1]);