function [] = PlotQ8Nodal(X, C, U)


[s,t]=meshgrid(-1:0.1:1,-1:0.1:1);

nElem = size(C, 1);


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
    surf(xx, yy, res ,'FaceColor', 'interp' , 'EdgeColor', 'none')
    hold on;
end

view(0,90)
colormap jet