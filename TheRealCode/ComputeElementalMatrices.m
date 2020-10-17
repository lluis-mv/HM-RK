% Function to initialize GPInfo and compute all matrices (B, N, H, ....)

function [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP)



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
    
    X = Nodes(Elements(el,:),:);
    % Linear triangles
    alfa = 1/3;
    beta = 1/3;
    
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
    Nsmall_chi = [-1 -1; 1 0; 0 1];
    Nu = (zeros(ndim, nNodes*ndim));
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
    
    GPInfo(el).Weight = Area;
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
    
end
