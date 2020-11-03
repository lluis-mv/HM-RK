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

xx = sym('x', [18,1], 'real')
Na = [a,a,a,b,b,b]';
Na = xx(1:6) + xx(7:12)*alfa+xx(13:18)*beta;

term = Na*(N-Na)'

I = int(term, beta, [0, 1-alfa]);
I = int(I, alfa, [0, 1])


I = reshape(I, 36,1);


aa = linspace(0,1,10);
bb  = linspace(0,1,10);
[aa,bb] = meshgrid(aa,bb);
for i = 1:size(aa,1)
    for j = 1:size(bb,1)
        if ( bb(i,j) < 1-aa(i,j))
            m = subs( sum(Na), alfa, aa(i,j));
            m = subs(m, beta,  bb(i,j));
            I = [I; (m-1)];
        end
    end
end

for i = 1:size(I,1)
    for j = 1:length(xx)
        J(i,j) = diff( I(i), xx(j));
    end
end

I = matlabFunction(I, 'file', 'residual');
J = matlabFunction(J, 'file', 'jacobian');

return;


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
    
