function [f, uDir, AllZero] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(3*nNodes, 1);


E = CP.E;
nu = CP.nu;

for el = 1:nElements
    Cel = Elements(el,:);
    
    for gp = 1:size(GPInfo, 2)
        
        Xpg = GPInfo(el, gp).Nu(1,1:2:end)*Nodes(Cel,:);
        y = Xpg(2);
        
        ff = (E*(nu - 1)*(6*y^2 - 6*y + 1))/(10*nu^2 + 5*nu - 5);
        ff2 = 0;

        ff = [0;ff];
        f( GPInfo(el, 1).dofsU) = f( GPInfo(el, 1).dofsU) - GPInfo(el, gp).Weight*GPInfo(el, gp).Nu'*( ff);
        
        f( GPInfo(el, 1).dofsWP) = f( GPInfo(el, 1).dofsWP) - GPInfo(el, gp).Weight*GPInfo(el, gp).N'*ff2;
    end
    
end


uDir = 0*f;
AllZero = false;