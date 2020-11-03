function NewIdea

syms alfa real
syms betta real

N =  [ (1 - alfa - betta)*(1-2*alfa-2*betta);
    alfa*(2*alfa-1);
    betta*(2*betta-1);
    4*alfa*(1-alfa-betta);
    4*alfa*betta;
    4*betta*(1-alfa-betta)];

a = [1, alfa, betta]';

De = IntegrateThis(a*a');
Ee = IntegrateThis(a*N');


Projection = Ee'*inv(De)*a;

this = IntegrateThis( Projection*(N-Projection)')

function RES = IntegrateThis(f)
syms betta real
syms alfa real
RES = int(f, betta, [0, 1-alfa]);
RES = int(RES, alfa, [0,1]);
