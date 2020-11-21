function IntegrateMassMatrix()


syms alfa real
syms betta real

N =  [ (1 - alfa - betta)*(1-2*alfa-2*betta);
    alfa*(2*alfa-1);
    betta*(2*betta-1);
    4*alfa*(1-alfa-betta);
    4*alfa*betta;
    4*betta*(1-alfa-betta)];

term = N*N';
term = [diff(N, alfa), diff(N, betta)];
term = term'*N;
Ma = IntegrateThis(term);



w = 1;
al = [1/3];
be = [1/3];

M1 = IntegrateNumThis(term, al, be, w);


w = 1/3*ones(1,3)*0.5;
al = [1/2, 1/2, 0];
be = [1/2, 0,1/3];
MZ3 = IntegrateNumThis(term, al, be, w);


w = 1/48*[-27, 25, 25, 25]*0.5;
al = [1/3, 0.6, 0.2, 0.2];
be = [1/3, 0.2, 0.6, 0.2];
MZ4 = IntegrateNumThis(term, al, be, w);




w = 1/6*ones(1,3);    
al = [2/3, 1/6,1/6];
al = [1/6, 1/6,2/3];
M2 = IntegrateNumThis(term, al, be, w);

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

al = auxK(:,1)';
be = auxK(:,2)';
w = auxK(:,3)';
 MK = IntegrateNumThis(term, al, be, w);


hola = 1;



function RES = IntegrateThis(f)
syms betta real
syms alfa real
RES = int(f, betta, [0, 1-alfa]);
RES = int(RES, alfa, [0,1]);
RES = eval(RES);


function RES = IntegrateNumThis(term, al, be, w)

syms alfa real
syms betta real
RES = zeros(size(term));

for i = 1:length(w)
    aux = subs(term, alfa, al(i));
    aux = subs(aux, betta, be(i));
    RES = RES + w(i)*aux;
end
RES = eval(RES);