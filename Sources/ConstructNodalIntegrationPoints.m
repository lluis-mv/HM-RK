

function [GPNodes] = ConstructNodalIntegrationPoints(CP, Nodes, Elements, GPElements)

De = zeros(6,6);

E = CP.E;
nu = CP.nu;

for i = 1:3
    for j = 1:3
        if ( i==j)
            De(i,j) = (1-nu);
        else
            De(i,j) = nu;
        end
    end
    De(i+3, i+3) = (1-2*nu)/2;
end

De = E/(1+nu)/(1-2*nu) * De;


D = De([1,2,4], [1,2,4]);

for nod = 1:size(Nodes,1)
    [candidates,trash] = find( Elements == nod);
    GPNodes(nod,1).NeigElement = candidates;
    GPNodes(nod,1).NeigNodes = unique(Elements(candidates,:));
    GPNodes(nod,1).Weight = sum( [GPElements( GPNodes(nod).NeigElement).Weight])/3;
    GPNodes(nod).D = D;
    GPNodes(nod).D6 = De;
    GPNodes(nod).StressNew = GPElements(1).StressNew;
    GPNodes(nod).StressPrev = GPElements(1).StressPrev;

    GPNodes(nod).StrainNew = zeros(6,1);
    GPNodes(nod).StrainPrev = zeros(6,1);

    GPNodes(nod).HistoryNew = GPElements(1).HistoryNew;
    GPNodes(nod).HistoryPrev = GPElements(1).HistoryPrev;

    GPNodes(nod).MCC = GPElements(1,1).MCC;
    GPNodes(nod).VonMises = GPElements(1,1).VonMises;
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

for nod = 1:size(Nodes,1)
    Qnod = zeros( 2*length(GPNodes(nod).NeigNodes), length(GPNodes(nod).NeigNodes) );
    mIdentity = [1;1;0];
    for el = GPNodes(nod).NeigElement'
        weight = GPElements(el).Weight;

        Celem = Elements(el,:);

        N = zeros(1, length(GPNodes(nod).NeigNodes) );
        for i = 1:3
            index = find( GPNodes(nod).NeigNodes == Celem(i));
            N(index) = 7/108;
            if ( GPNodes(nod).NeigNodes(index) == nod)
                N(index) = 11/54;
            end
        end
        
        
        Qnod = Qnod +  (GPNodes(nod).B'*mIdentity)* N*weight;

    end
    GPNodes(nod).Q = Qnod;
end

