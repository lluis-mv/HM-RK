syms alfa real


% N =  [ (1 - alfa - betta)*(1-2*alfa-2*betta);
%     alfa*(2*alfa-1);
%     betta*(2*betta-1);
%     4*alfa*(1-alfa-betta);
%     4*alfa*betta;
%     4*betta*(1-alfa-betta)];

N = [ (1-alfa)*(1-2*alfa);
    4*alfa*(1-alfa)
    alfa*(2*alfa-1)];

figure(1)
hold off
for i = 1:3
    ezplot(N(i), [0,1])
    hold on
end

Na = [N(1), 0, N(2), 0, N(3), 0;
    0, N(1), 0, N(2), 0, N(3)];

NN = int(Na, alfa, [0,1]);