% integrate consolidation equation
hold off
syms m intenger
syms XX real;
syms z real;
syms TT positive;
pw = 0;

aux = pi/2*(2*m+1);
     % 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT)
term = 2/aux * sin( aux * XX) * exp( - aux^2 * TT);
pw = pw+term ;

expr= int(pw, XX)
expr2 =  subs(expr, XX, z) - subs(expr, XX, 1)

simplify(expr2)
expr3 = simplify(expr2)


% some analytical solution
T = linspace(-4,4,30);
T = 10.^T;
d = ones(size(T));
for i = 1:length(T)
    for n = 0:20
        aux = 2*n+1;
        term = -8/pi^2 /(aux)^2 * exp( - aux^2*pi^2*T(i)/4);
        d(i) = d(i) + term;
    end
end
semilogx(T, d, 'k', 'linewidth', 3)        
hold on



xx = 0;
d = (1-xx)*ones(size(T));

for i = 1:length(T)
    expr5 = subs(expr3, z, xx);
    expr5 = subs(expr5, TT, T(i));
    for t = 0:20
        expr4 = subs(expr5, m, t);
        d(i) = d(i) + eval(expr4)
    end
end
semilogx(T, d, 'r*-')
hold on




xx = 0.3;
d = (1-xx)*ones(size(T));

for i = 1:length(T)
    expr5 = subs(expr3, z, xx);
    expr5 = subs(expr5, TT, T(i));
    for t = 0:20
        expr4 = subs(expr5, m, t);
        d(i) = d(i) + eval(expr4)
    end
end
semilogx(T, d, 'b*-')        
hold off
    