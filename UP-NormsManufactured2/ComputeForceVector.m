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
        
        ff =  (E*(nu - 1))/(5*(2*nu - 1)*(nu + 1));
        ff2 = 0;

        ff = [0;ff];
        f( GPInfo(el, 1).dofsU) = f( GPInfo(el, 1).dofsU) - GPInfo(el, gp).Weight*GPInfo(el, gp).Nu'*( ff);
        
        f( GPInfo(el, 1).dofsWP) = f( GPInfo(el, 1).dofsWP) - GPInfo(el, gp).Weight*GPInfo(el, gp).N'*ff2;
    end
    
end


uDir = 0*f;
AllZero = false;


nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
y = 1;
load = -(E*y*(nu - 1))/(5*(2*nu - 1)*(nu + 1));

if (size(Elements,2) == 3)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = false;
        for i = 1:3
            for j = i+1:3
                if ( any(Cel(i) == nodesTop))
                    if ( any(Cel(j) == nodesTop))
                        found = true;
                        ii = i;
                        jj = j;
                    end
                end
            end
        end
        if (found)
            nodi = Cel(ii);
            nodj = Cel(jj);
            XX = Nodes(nodi,:)-Nodes(nodj,:);
            
            normal = [XX(2), -XX(1)];
            normal = normal/norm(normal);
            fe = load * 0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
            
            index = [ 3*(nodi-1)+[1,2], 3*(nodj-1)+[1,2]];
            f(index) = f(index) + fe;
            
        end
    end
elseif (size(Elements,2) == 6)
    for el = 1:nElements
        Cel = Elements(el,:);
        found = 0;
        indexs = [];
        for i = 1:6
            if ( any(Cel(i) == nodesTop))
                found = found+1;
                indexs = [indexs, Cel(i)];
            end
        end
        if (found == 3)
            nodi = (indexs(1));
            nodj = (indexs(2));
            XX = Nodes(nodi,:)-Nodes(nodj,:);
            
            
            normal = [XX(2), -XX(1)];
            normal = normal/norm(normal);
            ff = [ 1/6,   0, 1/6,   0, 2/3,   0;
                    0, 1/6,   0, 1/6,   0, 2/3]';
            fe = load * ff*normal'*norm(XX);

            
            index = [];
            for i = 1:3
                index = [index, 3*(indexs(i)-1)+[1,2]];
            end
            f(index) = f(index) + fe;
        end
    end
    
end