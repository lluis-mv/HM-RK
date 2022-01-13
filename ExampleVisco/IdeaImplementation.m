

tau = sym('t', [2,2], 'real')
l = sym('ll', [2,2], 'real')


term1 = tau*l

for i = 1:4
    TT(i,:) = [tau(1,1); tau(2,2); tau(1,2); tau(2,1)];
end
LL = [l(1,1); l(2,2); l(1,2); l(2,1)];
M = sym(zeros(4,4));

M(1,1) = tau(1,1);
M(1,3) = tau(2,1);

M(2,2) = tau(2,2);
M(2,4) = tau(1,2);

M(3,1) = tau(1,2);
M(3,3) = tau(2,2);

M(4,2) = tau(2,1);
M(4,4) = tau(1,1);


term = l*tau
M*LL


M2 = sym(zeros(4,4));
M2(1,1) = tau(1,1);
M2(1,3) = tau(1,2);

M2(2,2) = tau(2,2);
M2(2,4) = tau(2,1);

M2(3,2) = tau(1,2);
M2(3,4) = tau(1,1);

M2(4,1) = tau(2,1);
M2(4,3) = tau(2,2);


res = l*tau+ tau*l'

this = M*LL+M2*LL
this1(1) = res(1,1);
this1(2,1) = res(2,2);
this1(3,1) = res(1,2);
this1(4,1) = res(2,1);
this-this1

syms t real
t = t;
Q = [cos(t), sin(t);
    -sin(t), cos(t)];
X = sym('XX',[2,1], 'real')

x = Q*[X];
v = diff(x, t);
for i = 1:2
    for j = 1:2
        F(i,j) = diff(x(i), X(j));
        L(i,j) = diff(v(i), X(j));
    end
end
l = L*inv(F);



    

