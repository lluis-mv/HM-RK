syms x real
syms y real

u= 0.25*(1-x^2-y^2)
u = x^4
grad = [diff(u, x), diff(u, y)];

then = diff(grad(1),x) + diff(grad(2), y)