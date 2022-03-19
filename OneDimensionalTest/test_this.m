

syms x real
syms y real

syms M real
syms x1 real
syms x2 real
syms x3 real

M = 10;
k = 1E-5;

u = x^2*(x-1)^2;
p = x*(x-1);
f1 = M*diff(diff(u,x), x) - diff(p, x);
f2 = diff(u,x)+k*diff( diff(p,x),x);


f1 = simplify(f1);
f2 = simplify(f2);

L = (x3-x1);
x2 = 0.5*(x1+x3);
N = [ 2*(x-x3)*(x-x2), 2*(x-x2)*(x-x1), -4*(x-x1)*(x-x3)]/L^2;

fTerm1 = int(f1*N, x, x1, x3);
fTerm1 = simplify(fTerm1);

Np = [(x3-x)/L, (x-x1)/L];

syms u1 real
syms u2 real
syms u3 real
syms pw1 real
syms pw2 real


normU = int(  ( u-N*[u1; u3; u2])^2, x, x1, x3);
normU = simplify(normU);
normP = int(  ( p-Np*[pw1; pw2])^2, x, x1, x3);
normP = simplify(normP);


fTerm2 = int(f2*Np, x, x1, x3);
fTerm2 = simplify(fTerm2);

matlabFunction(u, x1, x3, 'File','SolutionU');
matlabFunction(p, x1, x3, 'File','SolutionP');

matlabFunction(fTerm1', x1, x3, 'File', 'f1X', 'Vars', [x1, x3]);
matlabFunction(fTerm2', x1, x3, 'File', 'f2X', 'Vars', [x1, x3]);


matlabFunction(normU, x1, x3, 'File', 'ErrorNormU', 'Vars', [x1, x3, u1, u3, u2]);
matlabFunction(normP, x1, x3, 'File', 'ErrorNormP', 'Vars', [x1, x3, pw1, pw2]);
