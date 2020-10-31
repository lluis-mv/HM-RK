% Function to initialize GPInfo and compute all matrices (B, N, H, ....)

function [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType)

if ( all(ElementType == 'T3T3') )
    [GPInfo] = ComputeElementalMatricesT3T3(Nodes, Elements, CP);
elseif ( all(ElementType == 'T6T6') )
    [GPInfo] = ComputeElementalMatricesT6T6(Nodes, Elements, CP);
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



for el = 1:nElements
    
    Cel = Elements(el,:);
    X = Nodes(Cel,:);
    
    
    dofs = [];
    dofsU = [];
    for ii = 1:length(Cel)
        dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
        dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
    end
    
    
    
    
    % Linear triangles
    alfa = 1/3;
    beta = 1/3;
    
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
        b = [-dN_dX(1,i), 0; 0, -dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
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
    
    GPInfo(el,1).Weight = Area;
    GPInfo(el).B =B;
    GPInfo(el).dN_dX = dN_dX;
    GPInfo(el).N = Nsmall';
    GPInfo(el).Nu = Nu;
    GPInfo(el).he = he;
    GPInfo(el).Ms = Ms;
    GPInfo(el).D = D;
    GPInfo(el).D6 = De;
    
    
    
    GPInfo(el).StressNew = zeros(6,1);
    GPInfo(el).StressPrev = zeros(6,1);
    
    GPInfo(el).StrainNew = zeros(6,1);
    GPInfo(el).StrainPrev = zeros(6,1);
    
    GPInfo(el).HistoryNew = 0;
    GPInfo(el).HistoryPrev = 0;
    
    GPInfo(el).MCC = false;
    
    GPInfo(el).dofs = dofs;
    GPInfo(el).dofsU = dofsU;
    GPInfo(el).IndexReorder = [1,2,7,3,4,8,5,6,9];
    
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



for el = 1:nElements
    for gp = 1:3
        
        Cel = Elements(el,:);
        X = Nodes(Cel,:);
        
        
        
        
        dofs = [];
        dofsU = [];
        for ii = 1:length(Cel)
            dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
            dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
        end
        
        
        if ( gp == 1)
            alfa = 2/3; beta = 1/6;
        elseif ( gp == 2)
            alfa = 1/6; beta = 1/6;
        elseif ( gp == 3)
            alfa = 1/6; beta = 2/3;
        end
       
        
        Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*beta-1);
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
            b = [-dN_dX(1,i), 0; 0, -dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
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
        
        Ms = zeros(6,6);
        
        GPInfo(el,gp).Weight = Area/3;
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
        GPInfo(el,gp).IndexReorder = [1,2,13, 3,4,14, 5,6,15, 7,8,16, 9,10,17, 11,12,18];
    end
end




function [GPInfo] = ComputeElementalMatricesT6T7(Nodes, Elements, CP)




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



for el = 1:nElements
    for gp = 1:4
        
        Cel = Elements(el,:);
        X = Nodes(Cel,:);
        
        
        
        
        dofs = [];
        dofsU = [];
        for ii = 1:length(Cel)
            dofs = [ dofs, (Cel(ii)-1)*3 + [1,2,3] ];
            dofsU = [ dofsU, (Cel(ii)-1)*3 + [1,2] ];
        end
        
        
        if ( gp == 1)
            alfa = 1/3; beta = 1/3;
        elseif ( gp == 2)
            alfa = 3/5; beta = 1/5;
        elseif ( gp == 3)
            alfa = 1/5; beta = 3/5;
        elseif ( gp == 4)
            alfa = 1/5; beta = 1/5;
        end
       
        
        Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*beta-1);
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
            b = [-dN_dX(1,i), 0; 0, -dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
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
        
        Ms = zeros(6,6);
        
        GPInfo(el,gp).Weight = Area*2;
        if ( gp == 1)
            GPInfo(el,gp).Weight = GPInfo(el,gp).Weight*(-9/32);
        else
            GPInfo(el,gp).Weight = GPInfo(el,gp).Weight*(25/96);
        end
            
            
        GPInfo(el,gp).B = B;
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
        GPInfo(el,gp).IndexReorder = [1,2,13, 3,4,14, 5,6,15, 7,8,16, 9,10,17, 11,12,18];
    end
end


