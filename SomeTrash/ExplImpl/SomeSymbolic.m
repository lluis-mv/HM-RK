function [] = SomeSymbolic()




[a, b, c] = GetRungeKutta(4)


syms y real
y = sym(1)
syms dt real


for i = 1:length(b)
    yS = y;
    for j = 1:i-1
        yS = yS + dt*a(i,j)*k(1,j);
    end
    k(1,i) = sourceTerm(yS);
end
yE = y
for i = 1:length(b)
    yE = yE+dt*b(i)*k(1,i);
end

yE = simplify(yE)



syms k1 real
syms k2 real
syms k3 real

[a,bk] = GetImplRungeKutta(-3)



k1a = solve( k1 - sourceTerm(y+dt*a(1,1)*k1 + dt*a(1,2)*k2 + dt*a(1,3)*k3), k1)
k2a = solve( k2 - sourceTerm(y+dt*a(2,1)*k1a + dt*a(2,2)*k2 + dt*a(2,3)*k3),  k2);
k1a = subs(k1a, k2, k2a);
k3a = solve( k3 - sourceTerm(y+dt*a(3,1)*k1a + dt*a(3,2)*k2a + dt*a(3,3)*k3),  k3);

k1a = subs(k1a, k3, k3a);
k2a = subs(k2a, k3, k3a);
y2 = y + dt*b(1)*k1a + dt*b(2)*k2a+dt*b(3)*k3a;

hola = 1;
function [dy] = sourceTerm(y)
dy = y;

function [a, b] = GetImplRungeKutta(RK)

if ( RK == -1)
    a = 1;
    b = 1;
elseif ( RK == -2)
    a = [0,0;
        1/2, 1/2];
    b = [1/2, 1/2];
elseif (RK == -3)
    a = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
        5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
        5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
    b = [5/18, 4/9, 5/18];
end