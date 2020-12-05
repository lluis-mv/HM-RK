% Ensamble elemtal matrices to create C and K
% Compute the stabilization factor

function [C] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM)

if (nargin == 8)
    AlphaStabM = 1;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = sparse(nSystem, nSystem);


one = [1,1,0]';


Idev = eye(6)-1/3*[1,1,1,0,0,0]'*[1,1,1,0,0,0];

for el = 1:nElements
    
    for ngp = 1:size(GPInfo,2)
        
        
        Ddev = Idev*GPInfo(el,ngp).D6;
        kke = GPInfo(el,ngp).B'*Ddev([1,2,4],[1,2,4])*GPInfo(el,ngp).B;
        
        Q = GPInfo(el,ngp).B'*one * GPInfo(el,ngp).N;
        
        
        Qt = GPInfo(el,ngp).N' * [1,1,1] * GPInfo(el,ngp).D6(1:3,[1:2,4])*GPInfo(el,ngp).B;
        H = -GPInfo(el,ngp).N'*GPInfo(el,ngp).N*3;
    
        
        he = sqrt( sum([GPInfo(el,:).Weight]));
   
        if (all(ElementType == 'T3T3'))
            AlphaStab = 0;
        elseif ( all(ElementType == 'T6T6'))
            AlphaStab = 0;
        elseif ( all(ElementType == 'T6T3'))
            AlphaStab = 0;
        else
            disp(ElementType)
            error('this element does not exist. yet')
        end
        

        
        AlphaStab = AlphaStab*AlphaStabM;
        AlphaStab = -800;

        Ms = GPInfo(el,ngp).Ms * AlphaStab;
        
        if ( all(ElementType=='T6T3') )
            Q2 = [Qt; zeros(3,12)];
            Ms = [Ms+H, zeros(3,3);
                -1,-1,0,2,0,0;
                 0, -1, -1, 0 ,2, 0;
                 -1, 0, -1, 0, 0, 2];
            Ce = [kke, Q, zeros(12,3); Q2, Ms];
            
        else
        
            Ce = [kke, Q; Qt, Ms+H];
            
        end
        
        aux = GPInfo(el,ngp).IndexReorder;

        
        Ce = Ce(aux,aux);

        index = GPInfo(el,ngp).dofs;
        C(index,index) =  C(index,index) + Ce*GPInfo(el,ngp).Weight;
    end
end

