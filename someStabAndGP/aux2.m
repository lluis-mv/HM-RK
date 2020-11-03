function [] = aux2()

global RES
RES = [];



x = [121*ones(3,1); zeros(15,1)];
% x = rand(18,1)
load('thisVar.mat')
yy = ComputeThisError(x)
fun = @(x) ComputeThisError(x);


y = ComputeThisError(x)

for i = 1:30
    [x, f] = fminsearch( fun, x);
%     x(7:18) = 0;
    f
    hola = x*6;
    x(7:18) = 0;
    
end
x
save('thisVar.mat', 'x')

function Y = ComputeThisError(xx)
xx(7:18) = 0;

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





A = zeros(6,6);
a = xx(1:6); b= xx(7:12); c = xx(13:18);



for i = 1:length(wW)
    alfa = wA(i);
    beta = wB(i);
    

        N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*alfa-1);
            beta*(2*beta-1);
            4*alfa*(1-alfa-beta);
            4*alfa*beta;
            4*beta*(1-alfa-beta)];
         Na = (a + b*alfa + c*beta);

         term = Na*(N-Na)';
         
         A = A + term*wW(i);
end
% A = A(4:6,:);
Y = norm(A);
% return;


x = linspace(0,1,3);
y = linspace(0,1,3);
[x,y] = meshgrid(x,y);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if ( y(i,j) < 1-x(i,j))
            this = EvalPol(xx, x(i,j), y(i,j))-1;
            Y = Y + norm(this);
        end
    end
end


function n = EvalPol(xx, alfa, betta)
aa = xx(1:6);
bb = xx(7:12);
cc = xx(13:18);
Na = aa+bb*alfa+cc*betta;
n = sum(Na);