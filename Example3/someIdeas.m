
syms q real
syms B real
syms d real
syms H real
syms dt real
syms df real
syms sigma0 real

X0 = sym('x0', [2,1], 'real')

f = [df; 0];
k = B*d*B

C = [k, q; q, 0];
K = [0,0; 0, H];

X = X0 + dt*inv(C)*(f-K*X0);
mm = inv(C)*(f-K*X0);
sigma = sigma0 + dt*d*B*mm(1)


nu = 3;
np = 3;

u0 = sym('uu', [nu,1], 'real')
p0 = sym('pw', [np,1], 'real')
X0 = 0*[u0; p0];
sigma0 = 0*sym( 'ss', [2, 1], 'real')
B = sym('b', [2, nu], 'real')

H = sym('h', [np,np], 'real')

d = sym('dd', [2,2], 'real')
f = sym('ff', [nu+np,1], 'real')
f(1) = 0;
k = B'*d*B;

q = sym('qq', [np,nu], 'real')

C = [k, q; q', 0*H];
K = 0*C;
K(nu+1:end, nu+1:end) = H;

% BC
dofs = [1,6];
C(dofs,:) = 0; 
C(dofs,dofs) = eye(size(dofs,2));
K(dofs,:) = 0;
f(dofs) = 0;

X = X0 + C\(f-K*X0);
U = X(1:nu)
P = X(nu+1:end)
sigma = sigma0+ d*B*U;

res = B'*sigma + q*P -f(1:3)
res(1) = 0;
simplify(res)


