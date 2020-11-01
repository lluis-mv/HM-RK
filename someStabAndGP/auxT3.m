function x = aux2()

global RES
RES = [];


x = 1/3*ones(9,1);
% load('thisVar.mat')
x(4:9)= 0;

yy = ComputeThisError(x)
fun = @(x) ComputeThisError(x);



for i = 1:30
    [x, f] = fminsearch( fun, x);
    
    f
    hola = x*6;
    
end
x
save('thisVar.mat', 'x')

function Y = ComputeThisError(xx)


wA = [2/3,1/6,1/6];
wB = [1/6,1/6,2/3];

A = zeros(3,3);
% xx(4:9) = 0;
a = xx(1:3); b= xx(4:6); c = xx(7:9);


M = zeros(3,3);
Ms = zeros(3,3);
for i = 1:3
    alfa = wA(i);
    beta = wB(i);
    
    
    N =  [ (1 - alfa - beta);
        alfa;
        beta];
    
    Na = a + b*alfa+c*beta;
    
    
    term = Na*(N-Na)';
    A = A + term/3;
    
    M = M + N*N'/3*0.5;
    Ms = Ms + Na*Na'/3*0.5;
    
end

alfa = 1/3; beta = 1/3;

N2 =  [ (1 - alfa - beta);
    alfa;
    beta];





Y = 1000*norm(A);



x = linspace(0,1);
y = linspace(0,1);
[x,y] = meshgrid(x,y);
for i = 1:size(x,1)
    for j = 1:size(y,1)
        if ( y(i,j) < 1-x(i,j))
            this = EvalPol(xx, x(i,j), y(i,j))-1;
            Y = Y + norm(this)^2;
        end
    end
end


function n = EvalPol(xx, alfa, betta)
aa = xx(1:3);
bb = xx(4:6);
cc = xx(7:9);
Na = aa+bb*alfa+cc*betta;
n = sum(Na);