clear all

x = linspace(0, 1,100);
nNodes = length(x);


K = sparse(nNodes, nNodes);

E = 100;

for i = 1:nNodes-1
    h = x(i+1)-x(i);
    B = [-1/h, 1/h];
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) - B'*E*B*h;
end

K(1,:) = 0; K(1,1) = 1;
f = sparse(nNodes,1);
f(end) = 1;

u = -K\f;

figure(1)
plot(x, u);



for i = 1:nNodes-1
    xPG(i) = mean(x(i:i+1));
    Stress(i) = E*B*u(i:i+1);
    Strain(i) = B*u(i:i+1);
end

figure(2)
plot(xPG, Stress);


figure(3)
plot(xPG, Strain);