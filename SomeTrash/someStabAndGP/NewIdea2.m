function NewIdea2

syms alfa real
syms betta real

N =  [ (1 - alfa - betta);
    alfa;
    betta];

a = [1]';

De = IntegrateThis(a*a');
Ee = IntegrateThis(a*N')

Projection = Ee*inv(De)*a;
hola = 1;

function RES = IntegrateThis(f)
syms betta real
syms alfa real
RES = int(f, betta, [0, 1-alfa]);
RES = int(RES, alfa, [0,1]);
