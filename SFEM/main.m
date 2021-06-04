% First attemt to do something similar to SFEM

function [] = main()
figure(1); hold off
figure(2); hold off
figure(11); hold off
figure(12); hold off
addpath('../')
% 1. Define the problem



nDofs = 3;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.2;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

perme = 1E-5;


eSize =0.02;

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
Ce = zeros(nDofs*nNodes,nDofs*nNodes);
Ke = zeros(nDofs*nNodes,nDofs*nNodes);
mIdentity = [1,1,0]';
for el = 1:nElements
    
    dofsU = GPElements(el).dofsU;
    dofswP = GPElements(el).dofsWP;
    
    % Kuu
    Ce(dofsU,dofsU) = Ce(dofsU,dofsU) + GPElements(el).B'*De*GPElements(el).B*GPElements(el).Weight;
    
    % Kuwp
    Ce(dofsU, dofswP) = Ce(dofsU, dofswP) + GPElements(el).B'* mIdentity *GPElements(el).N * GPElements(el).Weight;
    
    % Kwpu
    Ce(dofswP, dofsU) = Ce(dofswP, dofsU) + GPElements(el).N'*mIdentity'* GPElements(el).B * GPElements(el).Weight;
    
    % Kwpwp
    Ke(dofswP, dofswP) = Ke(dofswP, dofswP)  - GPElements(el).dN_dX' * perme * GPElements(el).dN_dX * GPElements(el).Weight;
    
end


figure(1)
spy(Ke)

figure(2)
spy(Ce)



Cn = zeros(nDofs*nNodes, nDofs*nNodes);
Kn = zeros(nDofs*nNodes, nDofs*nNodes);
for nod = 1:nNodes
    CPatch = GPNodes(nod).NeigNodes;
    
    dofsU = [];
    dofswP = [];
    for mm = 1:length(CPatch)
        dofsU = [dofsU, (CPatch(mm)-1)*nDofs+[1, 2]];
        dofswP = [dofswP, (CPatch(mm)-1)*nDofs+3];
    end
    
    % Adding effective stresses Kuu
    Cn(dofsU,dofsU) = Cn(dofsU,dofsU) + GPNodes(nod).B'*De*GPNodes(nod).B*GPNodes(nod).Weight;
    
    Cn(dofsU, dofswP) =  Cn(dofsU, dofswP); 
    
    % Adding internal forces due to water phase Kuwp
    for el = GPNodes(nod).NeigElement'
        weight = GPElements(el).Weight*2;
        N = zeros(1, length(dofswP));
        Cel = Elements(el,:);
        for ind = 1:3
            index = find(dofswP == 3*(Cel(ind)-1)+3);
            N(index) = 7/216;
        end
        index = find(dofswP == 3*(nod-1)+3);
        N(index) = 11/108;
        
        Cn(dofsU, dofswP) = Cn(dofsU, dofswP) + GPNodes(nod).B'*mIdentity* N*weight;
    end
    
    % Adding mixsture deformation into mass balance Kwpu
    for el = GPNodes(nod).NeigElement'
        weight = GPElements(el).Weight*2;
        N = zeros(1, length(dofswP));
        Cel = Elements(el,:);
        for ind = 1:3
            index = find(dofswP == 3*(Cel(ind)-1)+3);
            N(index) = 7/216;
        end
        index = find(dofswP == 3*(nod-1)+3);
        N(index) = 11/108;
        Cn(dofswP, dofsU) = Cn(dofswP, dofsU) + N'*mIdentity'*GPNodes(nod).B*weight;
    end
    
    Kn(dofswP, dofswP) = Kn(dofswP, dofswP) - GPNodes(nod).dN_dX'*perme*GPNodes(nod).dN_dX*GPNodes(nod).Weight;
    
end


figure(11)
spy(Kn)
figure(12)
spy(Cn)

dt = 1;

ue = SolveProblem(Ke, Ce, dt, Nodes, Elements, nDofs);
un = SolveProblem(Kn, Cn, dt, Nodes, Elements, nDofs);
une = SolveProblem(Kn, Ce, dt, Nodes, Elements, nDofs);
uen = SolveProblem(Ke, Cn, dt, Nodes, Elements, nDofs);

figure(22)
plot(Nodes(:,2), ue(2:nDofs:end), 'k*')
hold on
plot(Nodes(:,2), un(2:nDofs:end), 'r*')
% plot(Nodes(:,2), une(2:nDofs:end), 'g*')
% plot(Nodes(:,2), uen(2:nDofs:end), 'b*')
legend('ue', 'un', 'une', 'uen', 'location', 'best')
hold off
legend('ue', 'un', 'une', 'uen')
hola = 1;
figure(23)
plot(Nodes(:,2), ue(3:nDofs:end), 'k*')
hold on
plot(Nodes(:,2), un(3:nDofs:end), 'r*')
% plot(Nodes(:,2), une(3:nDofs:end), 'g*')
% plot(Nodes(:,2), uen(3:nDofs:end), 'b*')
hold off
legend('ue', 'un', 'une', 'uen', 'location', 'best')
hola = 1;

figure(900); clf
subplot(1,2,1)
hold off;
PlotEffectiveStress(ue, Nodes, Elements, De, GPElements, nDofs)

subplot(1,2,2)
hold off;
PlotEffectiveStress(ue, Nodes, Elements, De, GPElements, nDofs, true)



figure(901); clf;
subplot(1,2,1)
hold off
PlotEffectiveStressNodal(un, Nodes, Elements, De, GPElements, GPNodes, nDofs);

subplot(1,2,2)
hold off
PlotEffectiveStressNodal(un, Nodes, Elements, De, GPElements, GPNodes, nDofs, true);







function [] = PlotEffectiveStressNodal(U, Nodes, Elements, De, GPElements, GPNodes, nDofs, total)

if ( nargin == 7)
    total = false;
end

nNodes = size(Nodes,1);
SigmaV = zeros(nNodes,1);

for nod = 1:nNodes
    NeigNodes = GPNodes(nod).NeigNodes;
    
    dofs = [];
    for ii = 1:length(NeigNodes)
        dofs = [dofs, nDofs*(NeigNodes(ii)-1)+1, nDofs*(NeigNodes(ii)-1)+2];
    end
    
    Upatch = U(dofs);
    Sigma = De*GPNodes(nod).B*Upatch;
    SigmaV(nod) = Sigma(2);
    
end

if ( total)
    SigmaV = SigmaV+U(3:3:end);
end

trisurf(Elements, Nodes(:,1), Nodes(:,2), SigmaV)

view(0,90)
axis equal
axis off
colorbar


function [] = PlotEffectiveStress(U, Nodes, Elements, De, GPElements, nDofs, total)

if (nargin == 6)
    total = false;
end


nElem = size(GPElements,1);
nNodes = size(Nodes,1);

XX = zeros(nElem,3);
YY = zeros(nElem,3);
EffVerticalStress = zeros(nElem,1);
WaterPressure = zeros(nElem,1);

for el = 1:nElem
    Ce = Elements(el,:);
    xE = Nodes(Ce,:);
    XX(el,:) = xE(:,1)';
    YY(el,:) = xE(:,2)';
    
    dofs = [];
    for ii = 1:3
        dofs = [dofs, nDofs*(Ce(ii)-1)+1, nDofs*(Ce(ii)-1)+2];
    end
    Uel = U(dofs);
    strain = GPElements(el).B*Uel;
    stress = De*strain;
    EffVerticalStress(el) = stress(2);
    dofswP = 3*(Ce-1)+3;
    WaterPressure(el) = 1/3*sum(U(dofswP));
end
if (total)
    EffVerticalStress = EffVerticalStress + WaterPressure; 
end


patch(XX',YY',EffVerticalStress);
axis equal
axis off
colorbar



function u = SolveProblem(K,C, dt, Nodes, Elements, nDofs)


nNodes = size(Nodes, 1);

f = ComputeForceVector(Nodes, Elements);

A = K + C/dt;


nodesBottom = find(Nodes(:,2) == 0);
nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));
nodesLeft = find(Nodes(:,1) == min(Nodes(:,1)));
nodesRight = find(Nodes(:,1) == max(Nodes(:,1)));




dofs = nDofs*(nodesBottom-1)+2;
A(dofs,:) = 0;
A(dofs,dofs) =  eye(length(dofs));

dofs = nDofs*([nodesBottom; nodesTop]-1)+1;
A(dofs,:) = 0;
A(dofs,dofs) =  eye(length(dofs));

dofs = nDofs*( [nodesTop]-1) + 3;
% dofs = nDofs*( [1:nNodes]'-1) + 3;
A(dofs,:) = 0;
A(dofs,dofs) =  eye(length(dofs));



dofs = nDofs*( [nodesLeft; nodesRight]-1) + 1;
A(dofs,:) = 0;
A(dofs,dofs) =  eye(length(dofs));

u = A\f;


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


for nod = 1:size(Nodes,1)
    dN_dXnod = zeros(2,  length(GPNodes(nod).NeigNodes));
    
    for neigElem = [GPNodes(nod).NeigElement]'
        dN_dXe = GPElements(neigElem).dN_dX  * GPElements(neigElem).Weight/3.0/GPNodes(nod).Weight;
        Celem = Elements(neigElem,:);
        ind = [];
        for jj = 1:3
            tt = find(GPNodes(nod).NeigNodes == Celem(jj));
            ind = [ind, tt];
        end
        dN_dXnod(:,ind) = dN_dXnod(:,ind) + dN_dXe;
        
    end
    GPNodes(nod).dN_dX = dN_dXnod;
end


function f = ComputeForceVector(Nodes, Elements)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

nodesTop = find(Nodes(:,2) == max(Nodes(:,2)));

f = zeros(3*nNodes, 1);

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
        fe = 0.5*[1,0;0,1;1,0;0,1]*normal'*norm(XX);
        
        index = [ 3*(nodi-1)+[1,2], 3*(nodj-1)+[1,2]];
        f(index) = f(index) - fe;
        
    end
end