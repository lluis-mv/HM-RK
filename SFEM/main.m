% First attemt to do something similar to SFEM



function [] = ExampleOne()
close all
addpath('../')
% 1. Define the problem






CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.49;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




eSize = 0.1;

model = createpde(1);

dx = 0.4; dy = 1;

R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');




Nodes = mesh.Nodes';
Elements = mesh.Elements';
nNodes = size(Nodes,1);
nElements = size(Elements,1);

figure(33)
triplot( Elements, Nodes(:, 1), Nodes(:, 2))

[GPElements] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');
[GPNodes] = ConstructNodalIntegrationPoints(Nodes, Elements, GPElements);

De = GPElements(1).D;
% ensamble elemental
Ke = zeros(2*nNodes,2*nNodes);
for el = 1:nElements
    Ce = Elements(el,:);
    entry = [];
    for mm = 1:3
        entry = [entry, Ce(mm)*2+[-1, 0]];
    end
    Ke(entry,entry) = Ke(entry,entry) + GPElements(el).B'*De*GPElements(el).B*GPElements(el).Weight;
end


figure(1)
spy(Ke)

Kn = zeros(2*nNodes, 2*nNodes);
for nod = 1:nNodes
    CPatch = GPNodes(nod).NeigNodes;
    
    entry = [];
    for mm = 1:length(CPatch)
        entry = [entry, CPatch(mm)*2+[-1, 0]];
    end
    
    Kn(entry,entry) = Kn(entry,entry) + GPNodes(nod).B'*De*GPNodes(nod).B*GPNodes(nod).Weight;
     
%       figure(1); clf; spy(Ke, 'r'); hold on;  pause(1); spy(Kn, 'b')
end


figure(2)
spy(Kn)


ue = SolveProblem(Ke, Nodes);
un = SolveProblem(Kn, Nodes);

figure(22)
plot(Nodes(:,2), ue(2:2:end), 'r*')
hold on
plot(Nodes(:,2), un(2:2:end), 'b*')

hola = 1;

function u = SolveProblem(K, Nodes)


nNodes = size(Nodes, 1);

f = zeros(2*nNodes, 1);





nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));


dofs = 2*(nodesTop-1)+2;
K(dofs,:) = 0;
f(dofs) = 0.001;
K(dofs,dofs) = eye(length(dofs));

dofs = 2*(nodesBottom-1)+2;
K(dofs,:) = 0;
K(dofs,dofs) =  eye(length(dofs));

dofs = 2*(nodesBottom-1)+1;
K(dofs,:) = 0;
K(dofs,dofs) =  eye(length(dofs));

u = K\f;


function [GPNodes] = ConstructNodalIntegrationPoints(Nodes, Elements, GPElements)

for nod = 1:size(Nodes,1)
    [candidates,trash] = find( Elements == nod);
    GPNodes(nod).NeigElement = candidates;
    GPNodes(nod).NeigNodes = unique(Elements(candidates,:));
    GPNodes(nod).Weight = sum( [GPElements( GPNodes(nod).NeigElement).Weight])/3;
end


for nod = 1:size(Nodes,1)
    Bnod = zeros(3,  2*length(GPNodes(nod).NeigNodes));
    
    for neigElem = [GPNodes(nod).NeigElement]'
        Be = GPElements(neigElem).B  * GPElements(neigElem).Weight/3.0/GPNodes(nod).Weight;
        Celem = Elements(neigElem,:);
        ind = [];
        for jj = 1:3
            tt = find(GPNodes(nod).NeigNodes == Celem(jj));
            ind = [ind, 2*tt-1, 2*tt];
        end
        Bnod(:,ind) = Bnod(:,ind) + Be;

    end
    
    GPNodes(nod).B = Bnod;
end

