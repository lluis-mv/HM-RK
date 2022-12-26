function [] = CreateTheTable()
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 2;


model = createpde(1);


R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';



ESIZE = 1./[2:11];



i = 1;

g = decsg(R1);
geometryFromEdges(model, g);

model1 = createpde(1);
geometryFromEdges(model1, g);




for eSize = ESIZE

    mesh = generateMesh(model, 'Hmax', eSize);
    Nodes = mesh.Nodes';
    Elements = mesh.Elements';



    mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');
    Nodes1 = mesh1.Nodes';
    Elements1 = mesh1.Elements';

    figure(1)
    triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'k')
    axis equal;
    axis off;
    print(['Mesh-', num2str(i)], '-dpdf');

    [GPInfo] = ComputeElementalMatrices(Nodes1, Elements1, CP, 'T3T3');
    he = mean(sqrt([GPInfo(:,:).Weight]));

    eSizeAxis(i) = he
    nElemens(i) = size(Elements1,1)
    nPatch(i) = size(Nodes1,1)
    nDofs(i) = size(Nodes1,1)*3
    nDofsQ(i) = size(Nodes,1)*3
    i = i+1;

end



eSize= 0.20;


figure(1); clf;
[Elements1, Nodes1] = CreateMesh(0,4, -4, 0, eSize, eSize);

triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'b')
[Nodes2] = AntiLaplacianSmoothing(Elements1, Nodes1);








for this = [-0.25,  0.25]



    Nodes1 = Nodes2;

    for j = 1:10
        [Elements1, Nodes1] = AntiLaplacianSmoothing(Elements1, Nodes1, this);
        figure(2); clf;
        triplot(Elements1, Nodes1(:,1), Nodes1(:,2), 'r')
    end

        [Nodes, Elements] = ConvertQuadratic(Nodes1, Elements1);
     [GPInfo] = ComputeElementalMatrices(Nodes1, Elements1, CP, 'T3T3');
    he = mean(sqrt([GPInfo(:,:).Weight]));

    eSizeAxis(i) = he
    nElemens(i) = size(Elements1,1)
    nPatch(i) = size(Nodes1,1)
    nDofs(i) = size(Nodes1,1)*3
    nDofsQ(i) = size(Nodes,1)*3
    i = i+1;
    
end

hola = 1;




function [C, X ] = AntiLaplacianSmoothing(C, X, this)
nNodes = size(X,1);
IsNeig = zeros(nNodes,1);

ind = find(X(:,1) == min(X(:,1)));
IsNeig(ind) = 1;

ind = find(X(:,1) == max(X(:,1)));
IsNeig(ind) = 1;
ind = find(X(:,2) == min(X(:,2)));
IsNeig(ind) = 1;
ind = find(X(:,2) == max(X(:,2)));
IsNeig(ind) = 1;
xL = X;
if (nargout == 1)
    for i = 1:nNodes
        if ( IsNeig(i) ==0)
            X(i,:) =  X(i,:) + 0.05*(rand(1,2)-0.5);
        end
    end
    C = X;
    return;
end

for i = 1:nNodes
    if ( IsNeig(i) ==0)
        [a,b] = find(C == i);
        %         NeigNodes = unique(C(a,:));
        NeigNodes = sort(C(a,:));
        [index] = find( NeigNodes~=i);
        NeigNodes = NeigNodes(index);
        xL(i,:) = mean(X(NeigNodes,:));
    end
end
displ = xL -X;
for i = 1:nNodes
    if ( norm(displ(i,:) )> 0)
        X(i,:) = X(i,:) +this *displ(i,:); %/norm(displ(i,:));
    end
end
C = delaunay(X(:,1), X(:,2));

function [C, X] = CreateMesh(xMin, xMax, yMin, yMax, dx, dy)

ix = floor((xMax-xMin)/dx);
jx = floor((yMax-yMin)/dy);
dx = (xMax-xMin)/ix;
dy = (yMax-yMin)/jx;

ii = 1;

for j = 1:jx+1
    for i = 1:ix+1
        X(ii,:) = [ xMin+dx*(i-1), yMin+dy*(j-1)];
        ii = ii+1;
    end
end

C = delaunay(X(:,1), X(:,2));
return;

nL = ix+1;
C = [];
for j = 1:jx
    for i = 1:ix
        t = (j-1)*(ix+1) + i;
        C = [C; t, t+1, t+1+nL; t, t+nL, t+nL];
    end
end



function [Nodes, Elements] = ConvertQuadratic(Nodes, Elements)

nElem = size(Elements,1);
nNodes = size(Nodes,1);

nThis = nNodes;
Elements = [Elements, zeros(nElem, 3)];
for el = 1:nElem
    Cel = Elements(el,1:3);
    Xel = Nodes(Cel,:);
    cNew = [0,0,0];
    for nn = 1:3
        if ( nn == 1)
            this = 1:2;
        elseif ( nn == 2)
            this = 2:3;
        elseif (nn == 3)
            this = [3,1];
        else
            clear this;
        end
        xNew = mean(Xel(this,:));
        index = find( Nodes(:,1) == xNew(1) & Nodes(:,2) == xNew(2));
        if ( isempty(index))
            Nodes = [Nodes; xNew];
            nThis = nThis+1;
            cNew(nn) = nThis;
        else
            cNew(nn) = index;
        end
    end
    Elements(el,4:6) = cNew;
end


