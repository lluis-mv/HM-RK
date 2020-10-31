% Ensamble elemtal matrices to create C and K
% Compute the stabilization factor

function [C, K] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, dt, implicit, AlphaStabM, a, b)

if (nargin == 5)
    implicit = false;
    AlphaStabM = 1;
elseif (nargin == 6)
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
    
    for ngp = 1:size(GPInfo,2)
        kke = GPInfo(el,ngp).B'*GPInfo(el,ngp).D*GPInfo(el,ngp).B;
        Q = GPInfo(el,ngp).B'*one * GPInfo(el,ngp).N;
        H = -GPInfo(el,ngp).dN_dX'*perme*GPInfo(el,ngp).dN_dX;
    
        he1 = GPInfo(el,ngp).Weight
        he = sqrt(GPInfo(el,ngp).Weight)
   
        AlphaStab = 8*perme*dt/he^2;
        if ( implicit)
            AlphaStab = -0.65*perme*dt/he^2;
    
            if (AlphaStab < 0)
                AlphaStab = 0;
            end
        end
        if (nargin == 9)
            AlphaStab = a/ConstModulus;%+b*perme*dt/he^2;
        end

        AlphaStab = AlphaStab*AlphaStabM;


        Ms = GPInfo(el,ngp).Ms * AlphaStab;
        Ce = [kke, Q; -Q', Ms];


        Ke = [0*kke, 0*Q; 0*Q', H];

        aux = GPInfo(el,ngp).IndexReorder;

        Ke = Ke(aux,aux);
        Ce = Ce(aux,aux);

        index = GPInfo(el).dofs;
        K(index,index) =  K(index,index) + Ke*GPInfo(el).Weight;
        C(index,index) =  C(index,index) + Ce*GPInfo(el).Weight;
    end
    
end

