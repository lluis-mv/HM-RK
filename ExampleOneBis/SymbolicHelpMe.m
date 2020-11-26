syms alfa real
syms beta real

ndim = 2;

X = sym('x', [3,2],  'real')




Nsmall =  [ 1 - alfa - beta; alfa;  beta];
Nsmall_chi = [-1 -1; 1 0; 0 1];
Nu = sym(zeros(ndim, 3*ndim));
for i = 1:3
    for dd = 1:2
        Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
    end
end
J = Nsmall_chi'*X;
dN_dX = inv(J)*Nsmall_chi';

B = [];
for i = 1:3
    b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
    B = [B, b];
end

Area = [1 1 1;
    X(1,1) X(2,1) X(3,1);
    X(1,2) X(2,2) X(3,2)];
Area = det(Area)/2;

he = 0;
for i = 1:3
    aux = 0;
    for j = 1:2
        aux = aux + dN_dX(j,i);
    end
    he = he + abs(aux);
end
he = sqrt(2)*he;
he = 4/he;

Ms = 1/18*[2,-1,-1;-1,2,-1;-1,-1,2];