function [] = testThisRKIdea()

DT = 10.^linspace(-4,0);
for RK = [(1:8)]
    i = 1;
    for dt = DT
        [err(i)] = IntegrateThis(dt, RK);
        i = i+1;
    end
    loglog(DT, err, '*-.')
    hold on
end
hold off


function err = IntegrateThis(dt, n)


xn = 1;



[a, b] = GetRungeKutta(n);
%     k = sym(zeros(1,n));
k = 0*b;

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






xa = xn*exp(dt);
ya = xa;
%     ya = taylor(xa, dt, 'order', n+1);

err = norm(ya-xnew);
%     err = simplify(err)
%     err = subs(err, xn, 1)
%     err = subs(err, dt, 1e-3)
%     err = eval(err)



function  y = source(x)
y = x;