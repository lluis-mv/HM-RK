
function [] = ExampleFour()
addpath('../Sources')

CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 20;



H = -1.5;



nX = 10;
nY = floor(nX*abs(H));

x = [];
y = [];
for i = 1:nX
    for j = 1:nY
        x = [x, (i-1)/(nX-1)];
        y = [y, (j-1)/(nY-1)*H];
    end
end
figure(212); clf;
figure(312); clf;
figure(412); clf;

for mesh = [0, 0.1]

    if (mesh == 0)
        XNAME = 'A';
        fign = 0;
        fign2 = 200;
        SPEC = 'b';
        SPEC2 = 'r';
    else
        XNAME = 'B';
        fign = 100;
        fign2 = 300;
        SPEC = 'c-.';
        SPEC2 = 'm-.';
    end

    tri = delaunay(x,y);
    figure(2); clf;
    triplot(tri, x, y)
    axis equal

    C = tri';
    X = [x; y]';
    for i = 1:100
        [C, X] = ApplyLaplacianSmoothing(C, X, mesh);
        axis equal
    end

    Nodes = X;
    Elements = C;

    figure(1)
    triplot(Elements, Nodes(:,1), Nodes(:,2), 'k')
    axis equal;
    axis off;
    print(['BiaxialMesh-N-', XNAME], '-dpdf');


    nSteps = 500;
    dt = 0.03/nSteps;


    tic
    [U, GPInfo, GPNodes, rrr,  information2] = ComputeImplicitNonLinearProblemNodal(Nodes, Elements, CP, dt, nSteps, 'T3T3', 1);
    toc
    FF = [information2.F];
    FF(1:2:end) = FF(1:2:end);


    figure(212);
    plot( [information2.t], FF(1:2:end), [SPEC], 'linewidth', 2,'DisplayName', ['NS-T3T3. ' XNAME])
    hold on
    

    figure(312);
    plot( [information2.t], FF(1:2:end), [SPEC], 'linewidth', 2, 'DisplayName', [XNAME])
    hold on
    

    figure(fign+21); clf
    SV = [GPNodes.StressNew];
    SV = SV(2,:);
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, SV);
    drawnow


    figure(fign+22); clf
    SV = [GPNodes.StressNew];
    pEff = mean(SV(1:3,:));
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, pEff);
    drawnow



    figure(fign+23); clf
    GPNodes = ComputeStrainInvatiants(GPNodes);
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, [GPNodes.StrainDev]');
    drawnow



    tic
    [U, GPInfo, rrr,  information2] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', 1);
    toc
    FF = [information2.F];
    FF(1:2:end) = FF(1:2:end);


    figure(212);
    plot( [information2.t], FF(1:2:end), [SPEC2], 'linewidth', 2,'DisplayName', ['T3T3. ' XNAME])
    hold on
    

    figure(412);
    plot( [information2.t], FF(1:2:end), [SPEC2], 'linewidth', 2, 'DisplayName', [XNAME])
    hold on
    

    figure(fign2+21); clf
    SV = [GPNodes.StressNew];
    SV = SV(2,:);
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, SV);
    drawnow


    figure(fign2+22); clf
    SV = [GPNodes.StressNew];
    pEff = mean(SV(1:3,:));
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, pEff);
    drawnow



    figure(fign2+23); clf
    GPNodes = ComputeStrainInvatiants(GPNodes);
    PlotHistoryVariableNodal( Nodes, Elements, GPNodes, [GPNodes.StrainDev]');
    drawnow





end














function [C, X ] = ApplyLaplacianSmoothing(C, X, this)
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





