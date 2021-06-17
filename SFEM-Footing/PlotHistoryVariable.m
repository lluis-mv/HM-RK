function PlotHistoryVariable( Nodes, C, GPNodes, values)


x = Nodes(:,1);
y = Nodes(:,2);



nElem = size(C,1);

xP = [];
yP = [];
c = [];

for elem = 1:nElem
        
    value = values(elem);
    
    Celem = C(elem,:);
    
    xP = [xP, x(Celem)];
    yP = [yP, y(Celem)];
    c = [c, value];
    

    
end
patch(xP, yP, c, 'edgecolor', 'none')

