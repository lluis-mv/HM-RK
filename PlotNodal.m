function [] = PlotNodal(X, C, U, showMesh)

if (nargin == 3)
    showMesh = false;
end

if ( size(C, 2) == 8)
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
elseif ( size(C, 2) == 6)
    [s,t]=meshgrid(0:0.05:1,0:0.05:1);
    for ii = 1:size(t,2)
        t(:,ii) = t(:,ii) * (1- s(end,ii));
    end
    
    nElem = size(C, 1);
    
    for elem = 1:nElem
        res = 0*t;
        for i = 1:size(s,1)
            for j = 1:size(t)
                alfa = s(i,j);
                beta = t(i,j);
                N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
                    alfa*(2*alfa-1);
                    beta*(2*beta-1);
                    4*alfa*(1-alfa-beta);
                    4*alfa*beta;
                    4*beta*(1-alfa-beta)];
                xx(i,j) = N'*X(C(elem,:), 1);
                yy(i,j) = N'*X(C(elem,:), 2);
                res(i,j) = N'*U(C(elem,:));
                
            end
        end
%         [a] = find( t > 1-s);
%         res(a) = nan;
        surf(xx, yy, res ,'FaceColor', 'interp' , 'EdgeColor', 'none')
        hold on;
    end
    
end

view(0,90)
colormap jet


if ( showMesh)
    if ( size(C, 2) == 8)
        tt = [1, 5, 2, 6, 3, 7, 4, 8, 1];
    elseif ( size(C, 2) == 6)
        tt = [1, 4, 2, 5, 3, 6, 1];
    end
    for elem = 1:nElem
        Cel = C(elem,:);
        plot3( X( Cel(tt), 1), X( Cel(tt), 2), 1E9*ones( size(X( Cel(tt), 2))),  'k*-')
        hold on;
    end
end