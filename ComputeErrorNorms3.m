
function [L2, L2U, LInf, LInfU] = ComputeErrorNorms3(X, Xa, Nodes, Elements, GPInfo, t, CP)

if (any(isnan(X)) ||any(isnan(Xa)) )
    L2 = nan; L2U= nan; LInf = nan; LInfU=nan;
    return
end
if (any(isinf(X)) ||any(isinf(Xa)) )
    L2 = nan; L2U= nan; LInf = nan; LInfU=nan;
    return
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


indexWP = sort(unique([GPInfo.dofsWP]));
LInf = max( abs( X(indexWP)-Xa(indexWP)));
LInfU = 0;
for i = 1:nNodes
    ind = 3*(i-1)+[1,2];
    thisNorm = norm(Xa(ind)-X(ind));
    LInfU = max(LInfU, thisNorm);
end


if ( size(Elements,2) == 8)
    [al, be, w] = GetWeights(9);
else
    [al, be, w] = GetWeights(6);
end

L2 = 0;
L2U = 0;


for el = 1:nElements
    
    dofsU = GPInfo(el,1).dofsU;
    dofsWP = GPInfo(el,1).dofsWP;
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:length(w)
        [Nu, Np] = GetShapeFunctions( al(gp), be(gp), length(dofsU), length(dofsWP));
        NN = Nu(1,1:2:end);
        xgp = NN * Nodes(Elements(el,:),:);
        u = 0;
        [v, pw] = EvaluateSticklePastor( 1-xgp(2),  t, 1, CP.M*CP.k, CP.M);
        v = -v;
        L2U = L2U + wA*w(gp)* norm( Nu*X(dofsU)- [u;v])^2;
        L2 = L2 + wA*w(gp)*abs( Np * X(dofsWP)-pw )^2;
    end
    
end

L2 = sqrt(L2/sum([GPInfo.Weight]));
L2U = sqrt(L2U/sum([GPInfo.Weight]));


if (L2 < 1E-15)
    L2 = rand*1E-15;
end
if (L2U < 1E-15)
    L2U = rand*1E-15;
end
if (LInf < 1E-15)
    LInf = rand*1E-15;
end
if (LInfU < 1E-15)
    LInfU = rand*1E-15;
end



function [Nu, Np] = GetShapeFunctions( alfa, beta, nU, nP)

ndim = 2;
if (nU == 6)
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
 
elseif (nU == 12)
    Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];    
elseif ( nU == 16)
    Nsmall =  [ -1/4*(1-alfa)*(1-beta)*(1+alfa+beta);
        -1/4*(1+alfa)*(1-beta)*(1-alfa+beta);
        -1/4*(1+alfa)*(1+beta)*(1-alfa-beta);
        -1/4*(1-alfa)*(1+beta)*(1+alfa-beta);
        1/2*(1-alfa^2)*(1-beta);
        1/2*(1+alfa)*(1-beta^2);
        1/2*(1-alfa^2)*(1+beta);
        1/2*(1-alfa)*(1-beta^2)];
elseif ( nU == 8)
    Nsmall =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)];
end

Nu = (zeros(ndim, size(Nsmall,1)*ndim));
for i = 1:size(Nsmall,1)
    for dd = 1:2
        Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
    end
end


if ( nP == 3)
    Np =  [ 1 - alfa - beta; alfa;  beta]';
elseif ( nP == 6)
    Np =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)]';
elseif ( nP == 4)
    Np =  1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)]';
elseif (nP == 8)
    Np =  [ -1/4*(1-alfa)*(1-beta)*(1+alfa+beta);
        -1/4*(1+alfa)*(1-beta)*(1-alfa+beta);
        -1/4*(1+alfa)*(1+beta)*(1-alfa-beta);
        -1/4*(1-alfa)*(1+beta)*(1+alfa-beta);
        1/2*(1-alfa^2)*(1-beta);
        1/2*(1+alfa)*(1-beta^2);
        1/2*(1-alfa^2)*(1+beta);
        1/2*(1-alfa)*(1-beta^2)]';
end

