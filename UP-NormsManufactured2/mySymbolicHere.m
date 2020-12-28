clear all

syms y real
syms x real
syms t real
syms E real
syms nu real


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



u =  0.1*y^2*t

epsilon = [diff(0,x); diff(u,y); 0; 0; 0; 0];

sigma = De *epsilon;
f = simplify(diff(   diff(sigma(2), y), t))

p = simplify(sum(sigma(1:3)))/3

alter = -simplify(diff(sigma(2), t))