

nNodes = 5;

xx = linspace(-1,1,nNodes)

syms chi real;
N = sym(ones(nNodes,1));

for i = 1:nNodes
    for j = 1:nNodes
        if ( i ~= j)
            N(i) = N(i)*(chi-xx(j))
        end
    end
end

for i = 1:nNodes
    a = subs(N(i),chi, xx(i));
    N(i) = N(i)/a;
end

hold off
for i = 1:nNodes
    ezplot(N(i), -1,1)
    ylim([-0.3,1.3])
    pause()
    hold on
end
plot(xx, ones(nNodes,1), '*')
plot(xx, zeros(nNodes,1), '*')
% 
%    (2*chi*(chi - 1)*(chi - 1/2)*(chi + 1/2))/3
%     -(8*chi*(chi - 1)*(chi + 1)*(chi - 1/2))/3
%  4*(chi - 1)*(chi + 1)*(chi - 1/2)*(chi + 1/2)
%     -(8*chi*(chi - 1)*(chi + 1)*(chi + 1/2))/3
%    (2*chi*(chi + 1)*(chi - 1/2)*(chi + 1/2))/3