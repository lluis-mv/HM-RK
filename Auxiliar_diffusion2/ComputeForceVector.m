function f = ComputeForceVector(t, Nodes, Elements, GPInfo)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(3*nNodes, 1);




for el = 1:nElements
    Cel = Elements(el,:);
    
    for gp = 1:size(GPInfo, 2)
        
        
        
        Xpg = GPInfo(el, gp).Nu(1,1:2:end)*Nodes(Cel,:);
        y = Xpg(2);
        
        ff = (2*t*(y - 1)^2)/5 - (t*y)/50 - (t*(y - 1))/50 + (2*t*y^2)/5 + (4*t*y*(2*y - 2))/5;
        ff2 = (2*t*y*(y - 1)^2)/5 - t^2/50 + (t*y^2*(2*y - 2))/5;

        ff = [0;ff];
        f( GPInfo(el, 1).dofsU) = f( GPInfo(el, 1).dofsU) - GPInfo(el, gp).Weight*GPInfo(el, gp).Nu'*( ff);
        
        f( GPInfo(el, 1).dofsWP) = f( GPInfo(el, 1).dofsWP) - GPInfo(el, gp).Weight*GPInfo(el, gp).N'*ff2;
    end
    
end
