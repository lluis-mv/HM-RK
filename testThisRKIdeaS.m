function [] = testThisRKIdeaS()
syms xn real
syms dt real

n = 5;

[a, b] = GetRungeKutta(n);

k = sym(0*b);

for i = 1:length(b)
    xStep = xn;
    for j = 1:i-1
        xStep = xStep + dt*a(i,j)*k(j);
    end
    k(i) = source(xStep);
end

xnew = xn;
for i = 1:length(b)
    xnew = xnew + dt*b(i)*k(i);
end

ii = sym(eye(1));

% Get ks
k2 = 0*k;
for i = 1:length(b)
    k2(i) = GetK(i, a, dt);
end
    
A = ii;
for i = 1:length(b)
    A = A + dt*b(i)*k2(i);
end

xnew2 = A*xn;

xa = xn*exp(dt);
%ya = xa;
ya = taylor(xa, dt, 'order', n+1);

err = norm(ya-xnew);
    err = simplify(err)
    err = subs(err, xn, 1)
    err = subs(err, dt, 1e-3)
    err = eval(err)

hola = 1;

function  y = source(x)
y = 8*x;


function k = GetK(i, a, dt)
num = 8;
if ( i == 1)
    k = num;
    return;
else
    k = num;
    for j = 1:i-1
        k = k + num*dt*a(i,j)*GetK(j, a, dt);
    end
end
