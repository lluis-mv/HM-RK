function [] = trash()
z = linspace(0, 1, 100);


for t = 10.^[-9:1:0]
pw = ComputeThis(z, 1, 1, t);
plot(z, pw) 
hold on
end
hold off

function [Xa] = ComputeThis(z, M, k, t)


% Other analytical solution...
nNodes = length(z);
Xa = zeros(nNodes,1);

for nod = 1:nNodes
    xx = 1-z(nod);
    TT = M * t*k;
    pw = 0;
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw = pw + 2/aux * sin( aux * xx) * exp( - aux^2 * TT);
    end
    Xa(nod) = pw;
end