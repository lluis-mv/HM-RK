Nodes1
Elements1

U = Nodes(:,1)+3*Nodes(:,2);

nElem = size(Elements,1);

X = Nodes;
C = Elements;
figure(1); clf;
[s,t]=meshgrid(-1:0.1:1,-1:0.1:1)
for elem = 1:nElem
    res = 0*t;
    for i = 1:size(s,1)
        for j = 1:size(t)
            alfa = s(i,j);
            beta = t(i,j);
            N =  [ -1/4*(1-alfa)*(1-beta)*(1+alfa+beta);
            -1/4*(1+alfa)*(1-beta)*(1-alfa+beta);
            -1/4*(1+alfa)*(1+beta)*(1-alfa-beta);
            -1/4*(1-alfa)*(1+beta)*(1+alfa-beta);
            1/2*(1-alfa^2)*(1-beta);
            1/2*(1+alfa)*(1-beta^2);
            1/2*(1-alfa^2)*(1+beta);
            1/2*(1-alfa)*(1-beta^2);
            ];
            xx(i,j) = N'*X(C(elem,:), 1);
            yy(i,j) = N'*X(C(elem,:), 2);
            res(i,j) = N'*U(C(elem,:));
        end
    end     
    surf(xx, yy, res)
    hold on;
end

