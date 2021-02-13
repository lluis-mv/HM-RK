function [] = main()

addpath('../')

nRK = 8;


[a, b, c] = GetRungeKutta(nRK);

a = sym('aa', [nRK, nRK], 'real');
b = sym('bb', [nRK, 1], 'real');


syms x0;
syms dt
syms B
for i = 1:length(a)
    xi = x0;
    for j = 1:i-1
        xi = xi+dt*a(i,j)*k(j)
    end
    k(i) = B*xi
end

xEnd = x0;
for i = 1:length(a)
    xEnd = xEnd + dt*b(i)*k(i);
end

A1 = diff(xEnd, x0)


A = ComputeAmplification(a, b, B, dt)

difere = simplify(xEnd-A*x0)

hola = 1;

% second Implementation

for i = 1:length(a)
    kk = 1;
    for j = 1:i-1
        kk = kk + dt*a(i,j)*ki(j);
    end
    ki(i) = B * kk
end
m = 1;
for i = 1:length(a)
    m = m + dt * b(i)*ki(i);
end

hola = 1;

function A = ComputeAmplification(a, b, B, dt)


A = 1;

for i = 1:length(a)
    A = A + dt*b(i)*B*GetThisTerm(a, b, B, dt, i);
end

function c = GetThisTerm(a, b, B, dt, i)

if ( i == 1)
    c = 1;
else
    c = 0;
    for j = 1:i-1
        c = c+ GetSubTerm( a, b, B, dt, i, j);
    end
end



function s = GetSubTerm(a, b, B, dt, i, j)

if ( j == 1)
    s = B*dt*a(i,j) + 1;
else
    s = B*dt*a(i,j)*GetSubTerm( a, b, B, dt, j, j-1);
end
