% Function to initialize GPInfo and compute all matrices (B, N, H, ....)

function [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType)

if ( all(ElementType == 'T3T3') )
    [GPInfo] = ComputeElementalMatricesT3T3(Nodes, Elements, CP);
elseif ( all(ElementType == 'T6T6') )
    [GPInfo] = ComputeElementalMatricesT6T6(Nodes, Elements, CP);
elseif ( all(ElementType == 'T6T3') )
    [GPInfo] = ComputeElementalMatricesT6T3(Nodes, Elements, CP);
end

hola = 1;

function [GPInfo] = ComputeElementalMatricesT3T3(Nodes, Elements, CP)




nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
ndim = 2;

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

al = 1/3;
be = 1/3;
w = 1;


for el = 1:nElements
    for gp = [1:length(w)]
        Cel = Elements(el,:);
        X = Nodes(Cel,:);
        
        
        dofs = [];
        dofsU = [];
        dofsWP = [];
        for ii = 1:length(Cel)
            dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
            dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
            dofsWP = [ dofsWP, ((Cel(ii)-1)*3 +3)];
        end
        
        
        
        
        % Linear triangles
        alfa = al(gp);
        beta = be(gp);
        
        Nsmall =  [ 1 - alfa - beta; alfa;  beta];
        Nsmall_chi = [-1 -1; 1 0; 0 1];
        Nu = (zeros(ndim, 3*ndim));
        for i = 1:3
            for dd = 1:2
                Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
            end
        end
        J = Nsmall_chi'*X;
        dN_dX = inv(J)*Nsmall_chi';
        
        B = [];
        for i = 1:3
            b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
            B = [B, b];
        end
        
        Area = [1 1 1;
            X(1,1) X(2,1) X(3,1);
            X(1,2) X(2,2) X(3,2)];
        Area = det(Area)/2;
        
        he = 0;
        for i = 1:3
            aux = 0;
            for j = 1:2
                aux = aux + dN_dX(j,i);
            end
            he = he + abs(aux);
        end
        he = sqrt(2)*he;
        he = 4/he;
        
        Ms = 1/18*[2,-1,-1;-1,2,-1;-1,-1,2];
        
        GPInfo(el,gp).Weight = Area*w(gp);
        GPInfo(el,gp).B =B;
        GPInfo(el,gp).dN_dX = dN_dX;
        GPInfo(el,gp).N = Nsmall';
        GPInfo(el,gp).Nu = Nu;
        GPInfo(el,gp).he = he;
        GPInfo(el,gp).Ms = Ms;
        GPInfo(el,gp).D = D;
        GPInfo(el,gp).D6 = De;
        
        
        
        GPInfo(el,gp).StressNew = zeros(6,1);
        GPInfo(el,gp).StressPrev = zeros(6,1);
        
        GPInfo(el,gp).StrainNew = zeros(6,1);
        GPInfo(el,gp).StrainPrev = zeros(6,1);
        
        GPInfo(el,gp).HistoryNew = 0;
        GPInfo(el,gp).HistoryPrev = 0;
        
        GPInfo(el,gp).MCC = false;
        
        GPInfo(el,gp).dofs = dofs;
        GPInfo(el,gp).dofsU = dofsU;
        GPInfo(el,gp).dofsWP = dofsWP;
        GPInfo(el,gp).dofsWPreal = dofsWP;
        GPInfo(el,gp).IndexReorder = [1,2,7,3,4,8,5,6,9];
    end
end




function [GPInfo] = ComputeElementalMatricesT6T6(Nodes, Elements, CP)




nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
ndim = 2;

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




al = [1/6,1/6,2/3];
be = [2/3,1/6,1/6];
w = [1/3,1/3,1/3];


for el = 1:nElements
    for gp = 1:length(w)
        
        Cel = Elements(el,:);
        X = Nodes(Cel,:);
        
        
        
        
        dofs = [];
        dofsU = [];
        dofsWP = [];
        for ii = 1:length(Cel)
            dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
            dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
            dofsWP = [ dofsWP, ((Cel(ii)-1)*3 + 3) ];
        end
        
        
        
        alfa = al(gp);
        beta = be(gp);
        
        Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*alfa-1);
            beta*(2*beta-1);
            4*alfa*(1-alfa-beta);
            4*alfa*beta;
            4*beta*(1-alfa-beta)];
        
        Nsmall_chi = [ 4*alfa + 4*beta - 3, 4*alfa + 4*beta - 3;
            4*alfa - 1,                   0;
            0,          4*beta - 1;
            4 - 4*beta - 8*alfa,             -4*alfa;
            4*beta,              4*alfa;
            -4*beta, 4 - 8*beta - 4*alfa];
        
        Nu = (zeros(ndim, 6*ndim));
        for i = 1:6
            for dd = 1:2
                Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
            end
        end
        J = Nsmall_chi'*X;
        dN_dX = inv(J)*Nsmall_chi';
        
        B = [];
        for i = 1:6
            b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
            B = [B, b];
        end
        
        
        Area = det(J)/2;
        
        he = 0;
        for i = 1:3
            aux = 0;
            for j = 1:2
                aux = aux + dN_dX(j,i);
            end
            he = he + abs(aux);
        end
        he = sqrt(2)*he;
        he = 4/he;
        
        Ms = 1/150*[   1, 1/3, 1/3,   -1,  1/3,   -1;
            1/3,   1, 1/3,   -1,   -1,  1/3;
            1/3, 1/3,   1,  1/3,   -1,   -1;
            -1,  -1, 1/3,  7/3, -1/3, -1/3;
            1/3,  -1,  -1, -1/3,  7/3, -1/3;
            -1, 1/3,  -1, -1/3, -1/3,  7/3];
        
        Ms = 1/360*[  6, -1, -1,  0, -4,  0;
            -1,  6, -1,  0,  0, -4;
            -1, -1,  6, -4,  0,  0;
            0,  0, -4, 12, -4, -4;
            -4,  0,  0, -4, 12, -4;
            0, -4,  0, -4, -4, 12];
        
        GPInfo(el,gp).Weight = Area*w(i);
        GPInfo(el,gp).B =B;
        GPInfo(el,gp).dN_dX = dN_dX;
        GPInfo(el,gp).N = Nsmall';
        GPInfo(el,gp).Nu = Nu;
        GPInfo(el,gp).he = he;
        GPInfo(el,gp).Ms = Ms;
        GPInfo(el,gp).D = D;
        GPInfo(el,gp).D6 = De;
        
        
        
        GPInfo(el,gp).StressNew = zeros(6,1);
        GPInfo(el,gp).StressPrev = zeros(6,1);
        
        GPInfo(el,gp).StrainNew = zeros(6,1);
        GPInfo(el,gp).StrainPrev = zeros(6,1);
        
        GPInfo(el,gp).HistoryNew = 0;
        GPInfo(el,gp).HistoryPrev = 0;
        
        GPInfo(el,gp).MCC = false;
        
        GPInfo(el,gp).dofs = dofs;
        GPInfo(el,gp).dofsU = dofsU;
        GPInfo(el,gp).dofsWP = dofsWP;
        GPInfo(el,gp).dofsWPreal = dofsWP;
        GPInfo(el,gp).IndexReorder = [1,2,13, 3,4,14, 5,6,15, 7,8,16, 9,10,17, 11,12,18];
    end
end



function [GPInfo] = ComputeElementalMatricesT6T3(Nodes, Elements, CP)




nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
ndim = 2;

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


wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;

wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;

auxK = [Na2, Na2, wa;
    Na1, Na2, wa;
    Na2, Na1, wa;
    Nb2, Nb2, wb ;
    Nb1, Nb2, wb ;
    Nb2, Nb1, wb ];

al = auxK(:,1)';
be = auxK(:,2)';
w = auxK(:,3)'/sum(auxK(:,3)');

for el = 1:nElements
    for gp = 1:length(w)
        
        Cel = Elements(el,:);
        X = Nodes(Cel,:);
        
        
        
        
        dofs = [];
        dofsU = [];
        dofsWP = [];
        for ii = 1:length(Cel)
            dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
            dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
            dofsWP = [ dofsWP, ((Cel(ii)-1)*3 + 3) ];
        end
        
        
        alfa = al(gp);
        beta = be(gp);
        
        NsmallP =  [ 1 - alfa - beta; alfa;  beta];
        Nsmall_chiP = [-1 -1; 1 0; 0 1];
        
        Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*alfa-1);
            beta*(2*beta-1);
            4*alfa*(1-alfa-beta);
            4*alfa*beta;
            4*beta*(1-alfa-beta)];
        
        Nsmall_chi = [ 4*alfa + 4*beta - 3, 4*alfa + 4*beta - 3;
            4*alfa - 1,                   0;
            0,          4*beta - 1;
            4 - 4*beta - 8*alfa,             -4*alfa;
            4*beta,              4*alfa;
            -4*beta, 4 - 8*beta - 4*alfa];
        
        Nu = (zeros(ndim, 6*ndim));
        for i = 1:6
            for dd = 1:2
                Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
            end
        end
        J = Nsmall_chi'*X;
        dN_dX = J\Nsmall_chi';
        
        B = [];
        for i = 1:6
            b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
            B = [B, b];
        end
        
        
        Area = det(J)/2;
        
        
        
        Ms = 1/18*[2,-1,-1;-1,2,-1;-1,-1,2];
        
        
        GPInfo(el,gp).Weight = Area*w(gp);
        GPInfo(el,gp).B =B;
        GPInfo(el,gp).dN_dX = J\Nsmall_chiP';
        GPInfo(el,gp).N = NsmallP';
        GPInfo(el,gp).Nu = Nu;
        
        GPInfo(el,gp).Ms = Ms;
        GPInfo(el,gp).D = D;
        GPInfo(el,gp).D6 = De;
        
        
        
        GPInfo(el,gp).StressNew = zeros(6,1);
        GPInfo(el,gp).StressPrev = zeros(6,1);
        
        GPInfo(el,gp).StrainNew = zeros(6,1);
        GPInfo(el,gp).StrainPrev = zeros(6,1);
        
        GPInfo(el,gp).HistoryNew = 0;
        GPInfo(el,gp).HistoryPrev = 0;
        
        GPInfo(el,gp).MCC = false;
        
        GPInfo(el,gp).dofs = dofs;
        GPInfo(el,gp).dofsU = dofsU;
        GPInfo(el,gp).dofsWPreal = dofsWP;
        GPInfo(el,gp).dofsWP = dofsWP(1:3); % Those relevant to compute internal forces,... blahblah
        GPInfo(el,gp).IndexReorder = [1,2,13, 3,4,14, 5,6,15, 7,8,16, 9,10,17, 11,12,18];
    end
end




