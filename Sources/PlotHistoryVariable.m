function PlotHistoryVariable( Nodes, C, GPNodes, values)


x = Nodes(:,1);
y = Nodes(:,2);

nElem = size(C,1);

xP = [];
yP = [];
c = [];

if ( size(values,2) == 1 )

    for elem = 1:nElem

        value = values(elem);

        Celem = C(elem,:);

        xP = [xP, x(Celem)];
        yP = [yP, y(Celem)];
        c = [c, value];

    end
    patch(xP, yP, c, 'edgecolor', 'none')

elseif ( size(values,2) == 9)
    [al, be, ~]  = GetWeights(9);
    A = [];
    for i = 1:length(al)
        alfa = al(i);
        beta = be(i);
        Nsmall =  [ -1/4*(1-alfa)*(1-beta)*(1+alfa+beta);
            -1/4*(1+alfa)*(1-beta)*(1-alfa+beta);
            -1/4*(1+alfa)*(1+beta)*(1-alfa-beta);
            -1/4*(1-alfa)*(1+beta)*(1+alfa-beta);
            1/2*(1-alfa^2)*(1-beta);
            1/2*(1+alfa)*(1-beta^2);
            1/2*(1-alfa^2)*(1+beta);
            1/2*(1-alfa)*(1-beta^2);
            ];
        A = [    A, Nsmall];
    end
    A = A';
    B = eye(4,8);
    nn = 1:4;


    [s,t]=meshgrid(-1:1:1,-1:1:1);
    %     [s,t]=meshgrid([-1,1],[-1,1]);

    for elem = 1:nElem

        v = values(elem,:);

        v3 = A\v';
        v2 = B*v3;

        Celem = C(elem,:);


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
                xx(i,j) = N'*Nodes(Celem, 1);
                yy(i,j) = N'*Nodes(Celem, 2);
                res(i,j) = N'*v3;
            end
        end
        for iii = 1:2
            surf(xx, yy, res ,'FaceColor', 'interp' , 'EdgeColor', 'interp')
            hold on;
        end

    end

    view(0,90)
else

    if ( size(values,2) == 3)
        A = [1, 1/6, 2/3;
            1,  1/6, 1/6;
            1, 2/3, 1/6];
        B = [1, 0,0;
            1, 1, 0;
            1, 0, 1];
        nn = 1:3;
    elseif ( size(values,2) == 4)


        [al, be, ~]  = GetWeights(4);
        A = [];
        for i = 1:length(al)
            alfa = al(i);
            beta = be(i);
            NsmallP =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)];
            A = [    A, NsmallP];
        end
        B = eye(4);
        nn = 1:4;

    end
        for elem = 1:nElem

            v = values(elem,:);

            v3 = A\v';
            v2 = B*v3;

            Celem = C(elem,nn);

            patch( x(Celem), y(Celem), v2, 'edgecolor', 'none')
        end


    end