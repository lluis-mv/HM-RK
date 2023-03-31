

function [GPNodes] = ConstructNodalIntegrationPointsQuad(CP, Nodes, Elements, GPElements)

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

ne = size(Elements,2)

for nod = 1:size(Nodes,1)
    [candidates,trash] = find( Elements == nod);
    GPNodes(nod,1).NeigElement = candidates;
    GPNodes(nod,1).NeigNodes = unique(Elements(candidates,:));

    ThisWeight = 0;
    for i = candidates'
        ind = find( Elements(i,:) == nod);
        ThisWeight = ThisWeight + GPElements(i,1).AllWeights(ind);
    end

    GPNodes(nod,1).Weight = ThisWeight;
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
        Celem = Elements(neigElem,:);
        cas = find( Celem == nod);
        [al, be, w]  = GetWeights(4);
        if ( cas == 1)
            al = al/2-0.5;
            be = be/2-0.5;
        elseif (cas == 2)
            al = al/2+0.5;
            be = be/2-0.5;
        elseif (cas == 3)
            al = al/2+0.5;
            be = be/2+0.5;
        elseif (cas == 4)
            al = al/2-0.5;
            be = be/2+0.5;
        end
        

        Be = 0*GPElements(neigElem,1).B;
        X = Nodes(Celem,:);
        for ii = 1:length(al)
            alfa = al(ii);
            beta = be(ii);
            weight = w(ii)*GPElements(neigElem,1).AllWeights(cas);

            B = [];
            Nsmall_chi = [   beta/4 - 1/4,   alfa/4 - 1/4;
                1/4 - beta/4, - alfa/4 - 1/4;
                beta/4 + 1/4,   alfa/4 + 1/4;
                - beta/4 - 1/4,   1/4 - alfa/4];
            J = Nsmall_chi'*X;
            dN_dX = inv(J)*Nsmall_chi';
            for i = 1:4
                b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
                B = [B, b];
            end
            Be = Be+B*weight;

        end





%         Be = GPElements(neigElem,1).B  * GPElements(neigElem,1).Weight/3.0/GPNodes(nod).Weight;
        
        Celem = Elements(neigElem,:);
        ind = [];
        for jj = 1:ne
            tt = find(GPNodes(nod).NeigNodes == Celem(jj));
            ind = [ind, 2*tt-1, 2*tt];
        end
        Bnod(:,ind) = Bnod(:,ind) + Be;
        
    end
    GPNodes(nod).B = Bnod / GPNodes(nod).Weight;
end


for nod = 1:size(Nodes,1)

    Qnod = zeros( 2*length(GPNodes(nod).NeigNodes), length(GPNodes(nod).NeigNodes) );
    mIdentity = [1;1;0];

    for neigElem = [GPNodes(nod).NeigElement]'
        Celem = Elements(neigElem,:);
        cas = find( Celem == nod);
        [al, be, w]  = GetWeights(4);
        if ( cas == 1)
            al = al/2-0.5;
            be = be/2-0.5;
        elseif (cas == 2)
            al = al/2+0.5;
            be = be/2-0.5;
        elseif (cas == 3)
            al = al/2+0.5;
            be = be/2+0.5;
        elseif (cas == 4)
            al = al/2-0.5;
            be = be/2+0.5;
        end
        

        Qe = zeros(8,4);
        X = Nodes(Celem,:);
        for ii = 1:length(al)
            alfa = al(ii);
            beta = be(ii);
            weight = w(ii)*GPElements(neigElem,1).AllWeights(cas);

            B = [];
            Nsmall_chi = [   beta/4 - 1/4,   alfa/4 - 1/4;
                1/4 - beta/4, - alfa/4 - 1/4;
                beta/4 + 1/4,   alfa/4 + 1/4;
                - beta/4 - 1/4,   1/4 - alfa/4];
            NsmallP =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)];
            J = Nsmall_chi'*X;
            dN_dX = inv(J)*Nsmall_chi';
            for i = 1:4
                b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
                B = [B, b];
            end
            Qe = Qe+B'*mIdentity*NsmallP'*weight;

        end

        
        Celem = Elements(neigElem,:);
        ind = [];
        ind2 = [];
        for jj = 1:ne
            tt = find(GPNodes(nod).NeigNodes == Celem(jj));
            ind = [ind, 2*tt-1, 2*tt];
            ind2 = [ind2, tt];
        end
        Qnod(ind, ind2) = Qnod(ind,ind2) + Qe;
        
    end
    GPNodes(nod).Q = Qnod / GPNodes(nod).Weight;
end




for nod = 1:size(Nodes,1)
    dN_dXnod = zeros(2,  length(GPNodes(nod).NeigNodes));
    dN_dXnod(1,:)= GPNodes(nod).B(1,1:2:end);
    dN_dXnod(2,:)= GPNodes(nod).B(2,2:2:end);

    GPNodes(nod).dN_dX = dN_dXnod;
end


