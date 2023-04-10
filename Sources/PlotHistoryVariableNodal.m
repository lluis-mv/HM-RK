function PlotHistoryVariableNodal( Nodes, C, GPNodes, values)


x = Nodes(:,1);
y = Nodes(:,2);



nNodes = length(x);

xP = [];
yP = [];
c = [];

if ( size(C,2) == 3)
    for node = 1:nNodes
        value = values(node);

        for elem = [GPNodes(node).NeigElement]'
            Celem = C(elem,:);
            index = find(Celem == node);
            if ( index == 1)
                i1 = Celem(2); i2 = Celem(3);
            elseif ( index == 2)
                i1 = Celem(1); i2 = Celem(3);
            else
                i1 = Celem(1); i2 = Celem(2);
            end

            x1 = [x(node), mean( [x(i1), x(node)]), mean(x(Celem)), mean([x(i2), x(node)])];
            y1 = [y(node), mean( [y(i1), y(node)]), mean(y(Celem)), mean([y(i2), y(node)])];
            xP = [xP, x1'];
            yP = [yP, y1'];
            c = [c, value];

        end
    end
elseif ( size(C,2) == 4)

    for node = 1:nNodes
        value = values(node);
        for elem = [GPNodes(node).NeigElement]'
            Celem = C(elem,:);
            index = find(Celem == node);
            if ( index == 1)
                i1 = Celem(2); i2 = Celem(4);
            elseif ( index == 2)
                i1 = Celem(3); i2 = Celem(1);
            elseif ( index == 3)
                i1 = Celem(4); i2 = Celem(2);
            elseif ( index == 4)
                i1 = Celem(1); i2 = Celem(3);
            end

            x1 = [x(node), mean( [x(i1), x(node)]), mean(x(Celem)), mean([x(i2), x(node)])];
            y1 = [y(node), mean( [y(i1), y(node)]), mean(y(Celem)), mean([y(i2), y(node)])];
            xP = [xP, x1'];
            yP = [yP, y1'];
            c = [c, value];

        end
    end
end
patch(xP, yP, c, 'edgecolor', 'none')

