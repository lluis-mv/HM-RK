syms u real
syms v real
syms x real
syms y real

E = 100;
nu = 0.3;


De = zeros(6,6);


for i = 1:3
    for j = 1:3
        if ( i==j)
            De(i,j) = (1-nu);
        else
            De(i,j) = nu;
        end
    end
    De(i+3, i+3) = (1-2*nu)/2;
end

De = E/(1+nu)/(1-2*nu) * De;


D = De([1,2,4], [1,2,4]);

u = 0.01*x^3*y;
v = 0.01/4*x^4;

u = 0;
v = y^3*(y-1)^2*0.01;
X = [x; y];





epsilon = [diff(u,x), diff(v,y), diff(u,y)+diff(v,x)]'
sigma = D*epsilon;

f1 = diff(sigma(1), x)+diff(sigma(3), y)
f2 = diff(sigma(2), y)+diff(sigma(3), x)


pw = (x-0.2)*x*y*(y-1);


p = De(1:3,1:3)*[epsilon(1:2);0];
pw = sum(p)/3;


f = diff(diff( pw, x), x) + diff( diff(pw, y), y);
f = 0;


matlabFunction( [-f1;-f2], x, y , 'File', 'source')
matlabFunction( u,v,pw, x, y , 'File', 'solution')

matlabFunction( [f], x, y , 'File', 'sourceWP')