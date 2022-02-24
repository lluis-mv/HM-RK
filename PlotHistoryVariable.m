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
    
else
    
       
    A = [1, 1/6, 2/3;
        1,  1/6, 1/6;
        1, 2/3, 1/6];
    B = [1, 0,0;
        1, 1, 0;
        1, 0, 1];
    
    for elem = 1:nElem
        
        v = values(elem,:);
        
        v3 = A\v';
        v2 = B*v3;
        
        Celem = C(elem,1:3);
        
        xP = [xP, x(Celem)];
        yP = [yP, y(Celem)];
        c = [c, v2];
        
    end
    surf(xP, yP, c, 'edgecolor', 'none')
    
end