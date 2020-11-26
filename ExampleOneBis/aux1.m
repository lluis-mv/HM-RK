function [GPInfo] = aux1()

Nodes=[0.0, 0.0;
    1, 0.0;
    0, 1]
Elements = [1,2,3];
CP.E = 1000; CP.nu = 0.2;

[GP1] = ComputeElementalMatricesT3T3(Nodes, Elements, CP, 1);
[GP2] = ComputeElementalMatricesT3T3(Nodes, Elements, CP, 3);
[GP3] = ComputeElementalMatricesT3T3(Nodes, Elements, CP, 3);
[GP4] = ComputeElementalMatricesT3T3(Nodes, Elements, CP, 4);

one = [1,1,0]';
Q1 = GP1.B'*one*GP1.N*GP1.Weight;

Q2 = 0*Q1;
for i = 1:length(GP2)
    Q2 = Q2 + GP2(i).B'*one*GP2(i).N*GP2(i).Weight;
end

Q3 = 0*Q1;
for i = 1:length(GP3)
    Q3 = Q3 + GP3(i).B'*one*GP3(i).N*GP3(i).Weight;
end


Q4 = 0*Q1;
for i = 1:length(GP4)
    Q4 = Q4 + GP4(i).B'*one*GP4(i).N*GP4(i).Weight;
end

hola = 1;

function [GPInfo] = ComputeElementalMatricesT3T3(Nodes, Elements, CP, GP)




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

if ( GP == 1)
al = 1/3;
be = 1/3;
w = 1;
elseif ( GP==3)
al = [1/6,1/6,2/3];
be = [2/3,1/6,1/6];
w = [1/3,1/3,1/3];
elseif (GP==4)
    al = [1/3, 3/5,1/5,1/5];
    be = [1/3, 1/5, 1/5,3/5];
    w = (1/48)*[-27, 25, 25, 25];
elseif (GP==6)
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
end

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


