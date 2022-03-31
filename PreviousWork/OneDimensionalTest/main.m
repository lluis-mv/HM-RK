function [] = main()

close all; clear all; clc; clf; 
test_this
close all; clear all; clc; clf; 

h = [2, 5, 10, 20, 50, 100, 200, 500];
for i = 1:length(h)
    [normU(i), normP(i)] = EvaluateProblem(h(i))
    figure(3)
    loglog( h(1:i), normU, 'r-.o',  h(1:i), normP, 'b-.o');
    drawnow;
    if ( i > 1)
        slopeU(i) = log10(normU(i)/normU(i-1))/log10(h(i)/h(i-1))
        slopeP(i) = log10(normP(i)/normP(i-1))/log10(h(i)/h(i-1))
    end
end
hola = 1;

function [normU, normP] = EvaluateProblem(nNodes)



M = 10;
k = 1E-5;
% nNodes = 2;
x = linspace(0,1,nNodes);
xpw = x;
C = [1:nNodes-1; 2:nNodes; nNodes+1:nNodes+nNodes-1];

xNew = [];
for el = 1:nNodes-1
    Cel = C(1:2,el);
    xNew = [xNew, mean(x(Cel))];
end
x = [x, xNew];



syms xx real
L = x(2)-x(1);

N = [ 2*(xx-L)*(xx-0.5*L),2*(xx-0.5*L)*xx, -4*xx*(xx-L), ]/L^2;
gradN = diff(N, xx);


Np = [(L-xx)/L, xx/L];
gradNp = diff(Np, xx);

K = -int( gradN'*M*gradN, xx, 0, L);

Q = int( gradN'*Np, xx, 0, L);

H = -k*int( gradNp'*gradNp, xx, 0, L);


K = eval(K); Q = eval(Q); H = eval(H); 

% ensamble
nNodesOr = nNodes;
nNodes = length(x);
nEl = size(C,2);
nNodesWP = nNodesOr;
A = zeros(nNodes+nNodesWP, nNodes+nNodesWP);
f = zeros(nNodes+nNodesWP, 1);
for el = 1:nEl
    Cel = C(:,el);
    iU = C(:,el);
    iPW = C(1:2,el)+nNodes;
    
    A(iU,iU) = A(iU,iU) + K;
    A(iU,iPW) = A(iU,iPW) +Q;
    A(iPW,iU) = A(iPW,iU) +Q';
    A(iPW,iPW) = A(iPW,iPW)+H; 
    f(iU) = f(iU) + f1X( x(Cel(1)), x(Cel(2)));
    f(iPW) = f(iPW) + f2X( x(Cel(1)), x(Cel(2)));
end

indexFix = [1, nNodesWP, nNodes+1,nNodes+nNodesWP];
A(indexFix,:) = 0;
A(indexFix,indexFix) = eye(length(indexFix), length(indexFix));
f(indexFix) = 0;
f(end) = 0;

u = A\f;

[xOrd, index ] = sort(x);
figure(1); clf;
subplot(2,1,1)

plot(xOrd, u(index), 'b', xOrd, SolutionU(xOrd), 'r-.')
subplot(2,1,2)
plot(xpw, u(nNodes+1:end), 'b', xpw, SolutionP(xpw), 'r-.')

% now lets compute the error norms
normU = 0; 
normP = 0;
for el = 1:nEl
    Cel = C(:,el);
    iU = C(:,el);
    iPW = C(1:2,el)+nNodes;
    
    normU = normU + ErrorNormU( x(Cel(1)), x(Cel(2)), u(Cel(1)), u(Cel(2)), u(Cel(3))) ;
    normP = normP + ErrorNormP( x(Cel(1)), x(Cel(2)), u(Cel(1)+nNodes), u(Cel(2)+nNodes)) ;
end
normU = real(sqrt(normU))
normP = real(sqrt(normP))

hola = 1;