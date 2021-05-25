function [] = Poisson()

NNODES = 2.^[2:2:10];
i = 1
for nNodes = NNODES
    
    [errE(i), errN(i)] = computeThis(nNodes);
    
    
    figure(222)
    loglog(NNODES(1:i), errE, 'r*-.', NNODES(1:i), errN, 'b*-.')
    
    
    i = i+1;
end


function [errE, errN] = computeThis(nNodes)

figure(1)

x = linspace(0, 1, nNodes)';
% x = sort(rand(1,nNodes))'
x(1) = 0; x(end) = 1;

nNodes = length(x);
Ke = zeros(nNodes, nNodes);
f = zeros(nNodes, 1);

for i = 1:nNodes-1
    h = x(i+1)-x(i);
    index = [i,i+1];
    
    Ke(index,index) =  Ke(index,index) + [1, -1; -1, 1]/h;
    f(index) = f(index) + 2 * [0.5; 0.5]*h;
end
f(1) = 0; f(nNodes) = 0;
Ke(1,:) = 0; Ke(1,1) = 1;
Ke(nNodes,:) = 0; Ke(nNodes,nNodes) = 1;
ue = -Ke\f;

% xa = linspace(0, 1, 1000);
ua = x.*(x-1)
plot(x, ua, 'k')
hold on
plot(x, ue, 'r')



NodalStructure = ConstructNodalStructure(x);
Kn = zeros(nNodes, nNodes);
for nod = 1:nNodes
    index = NodalStructure(nod).NeigNod;
    Kn(index,index) =  Kn(index,index) + NodalStructure(nod).grad'* NodalStructure(nod).grad*NodalStructure(nod).h;
end
Kn(1,:) = 0; Kn(1,1) = 1;
Kn(nNodes,:) = 0; Kn(nNodes,nNodes) = 1;
un = -Kn\f

figure(1)

plot(x, un, 'm')
hold off


% PlotTheGradients(x, NodalStructure)

errE = max(abs(ue-ua))
errN = max(abs(un-ua))


function PlotTheGradients(x, NodalStructure)

thisShape = 4;

figure(33);



for p = 1:length(NodalStructure)
    ind = find(NodalStructure(p).NeigNod == thisShape);
    if ( p > 1)
        xini = 0.5*(x(p)+x(p-1));
    else
        xini = 0;
    end
    if ( p < length(NodalStructure) )
        xfin = 0.5*(x(p)+x(p+1));
    else
        xfin = 1;
    end
    if ( length(ind) == 1)
        plot([xini, xfin], NodalStructure(p).grad(ind)*[1,1], 'k')
        hold on
    else
        plot([xini, xfin], [0,0], 'k')
        hold on
    end
end


if ( thisShape > 1)
    plot([x(thisShape-1), x(thisShape)], [1,1]/(x(thisShape)-x(thisShape-1)), 'r')
end
if ( thisShape < length(NodalStructure) )
    plot([x(thisShape+1), x(thisShape)], [-1,-1]/(x(thisShape+1)-x(thisShape)), 'r')
end
xlim([0, 1])
hold on


hold off


function [NodalStructure] = ConstructNodalStructure(Nodes)
nNodes = length(Nodes);
C = [1:nNodes-1; 2:nNodes]';

h = Nodes(2:end)-Nodes(1:end-1);

for nod = 1:nNodes
    [NeigElements, trash] = find(C == nod);
    NodalStructure(nod).NeigElem = NeigElements;
    NodalStructure(nod).NeigNod = unique( C(NeigElements,:));
    NodalStructure(nod).h = sum( h(NeigElements))/2; 
end

for nod = 1:length(Nodes)
    gradNod = zeros(1, length(NodalStructure(nod).NeigNod));
    
    for neigElem = [NodalStructure(nod).NeigElem]'
        grad = [-1,1]./h(neigElem);
        Celem = C(neigElem,:);
        ind  = grad;
        for jj = 1:2
            ind(jj) = find( NodalStructure(nod).NeigNod == Celem(jj));
        end
        gradNod(ind) = gradNod(ind) + grad* h(neigElem)/2/NodalStructure(nod).h;
    end
    NodalStructure(nod).grad = gradNod;
end