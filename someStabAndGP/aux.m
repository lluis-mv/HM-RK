
syms alfa real
syms beta real
syms h real

N = [1-alfa-beta, alfa, beta]';

Na = [1/3,1/3,1/3]'

term = N*N'-Na*Na';
term = Na*(N'-Na')

A = int(term, beta, [0, 1-alfa]);
A = int(A, alfa, [0, 1])

        N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*alfa-1);
            beta*(2*beta-1);
            4*alfa*(1-alfa-beta);
            4*alfa*beta;
            4*beta*(1-alfa-beta)];
        
%         aa = sym('a', [6,1], 'real')
%         bb = sym('b', [6,1], 'real')
%         cc = sym('c', [6,1], 'real')
xx = sym('x', [18,1], 'real')
syms a real
syms b real


        Na = xx(1:6)+xx(7:12)*alfa+xx(13:18)*beta; 
        Na = [a,a,a,b,b,b]';
%         Na = [1-alfa-beta, alfa, beta, 1/6,1/6,1/6]';
%         Na = 1/6*ones(6,1);

        
        term = Na*(N-Na)'


A = int(term, beta, [0, 1-alfa]);
A = int(A, alfa, [0, 1])
        
return

term2 = subs(term, alfa, 1/3);
term2 = subs(term2, beta, 1/3)/2

% now 1d

N = [1-alfa/h, alfa/h]';
Na = [1/2,1/2]';

term = N*N'-Na*Na'
A = int(term, alfa, [0, h])