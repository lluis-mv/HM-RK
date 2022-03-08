function [] = PlotMesh(X, C, plotNodes)

if ( nargin == 2)
    plotNodes = false;
end

nElem = size(C, 1);
if ( size(C, 2) == 8)
    tt = [1, 5, 2, 6, 3, 7, 4, 8, 1];
elseif ( size(C, 2) == 6)
    tt = [1, 4, 2, 5, 3, 6, 1];
end
for elem = 1:nElem
    Cel = C(elem,:);
    if (plotNodes )
        SPEC = 'k*-';
    else
        SPEC = 'k';
    end
    plot3( X( Cel(tt), 1), X( Cel(tt), 2), 1E9*ones( size(X( Cel(tt), 2))),  SPEC)
    hold on;
end
view(0,90)