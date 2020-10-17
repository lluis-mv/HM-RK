% Ensamble elemtal matrices to create C and K
% Compute the stabilization factor

function [C, K] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, AlphaStabM)

if (nargin == 5)
    AlphaStabM = 1;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = zeros(nSystem, nSystem);
K = C;


one = [1,1,0]';

perme = CP.k;
ConstModulus=CP.M;

for el = 1:nElements
    ind = Elements(el,:);
    index = [];
    for ii = 1:length(ind)
        index = [ index, (ind(ii)-1)*3 + [1,2,3] ];
    end
    
    kke = GPInfo(el).B'*GPInfo(el).D*GPInfo(el).B;
    Q = GPInfo(el).B'*one * GPInfo(el).N;
    H = -GPInfo(el).dN_dX'*perme*GPInfo(el).dN_dX;
    
    
    he = sqrt(GPInfo(el).Weight);
    %AlphaStab = 3/ConstModulus+6*perme*dt/he^2;
    
    %     AlphaStab = AlphaStab*2;
    
    AlphaStab = 8*perme*dt/he^2;
    AlphaStab = AlphaStab*AlphaStabM;
    
    
    Ms = GPInfo(el).Ms * AlphaStab;
    Ce = [kke, Q; -Q', Ms];
    Ke = [zeros(6,9); zeros(3,6), H];
    
    aux = [1,2,7,3,4,8,5,6,9];
    
    Ke = Ke(aux,aux);
    Ce = Ce(aux,aux);
    
    K(index,index) =  K(index,index) + Ke*GPInfo(el).Weight;
    C(index,index) =  C(index,index) + Ce*GPInfo(el).Weight;
    
end

