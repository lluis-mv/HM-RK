function [f, uDir, AllZero, F] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));

f = zeros(3*nNodes, 1);



load = 10;
load2 = 10*t;
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
            fe = load*0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
            
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
            fe = load*ff*normal'*norm(XX);

            
            index = [];
            for i = 1:3
                index = [index, 3*(indexs(i)-1)+[1,2]];
            end
            f(index) = f(index) + fe;
        end
    end
    
end

uDir = 0*f;
AllZero = true;

if (nargout < 4)
    return;
end




load = load2;
F = 0*f;
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
            fe = load*0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
            
            index = [ 3*(nodi-1)+[1,2], 3*(nodj-1)+[1,2]];
            F(index) = F(index) + fe;
            
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
            fe = load*ff*normal'*norm(XX);

            
            index = [];
            for i = 1:3
                index = [index, 3*(indexs(i)-1)+[1,2]];
            end
            F(index) = F(index) + fe;
        end
    end
    
end