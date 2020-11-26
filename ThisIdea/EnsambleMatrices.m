% Ensamble elemtal matrices to create C and K
% Compute the stabilization factor

function [C, K] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM, a, b)

if (nargin == 8)
    AlphaStabM = 1;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = sparse(nSystem, nSystem);
K = C;


one = [1,1,0]';

perme = CP.k;
ConstModulus=CP.M;

for el = 1:nElements
    
    for ngp = 1:size(GPInfo,2)
        kke = GPInfo(el,ngp).B'*GPInfo(el,ngp).D*GPInfo(el,ngp).B;
        Q = -GPInfo(el,ngp).B'*one * GPInfo(el,ngp).N;
        H = GPInfo(el,ngp).dN_dX'*perme*GPInfo(el,ngp).dN_dX;
    
        
        he = sqrt( sum([GPInfo(el,:).Weight]));
   
        if (all(ElementType == 'T3T3'))
        
            AlphaStab = 8*perme*dt/he^2;
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (8E4*dt*perme/he^2)^(2*RKMethod) ));
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (8E5*dt*perme/he^2)^(RKMethod) ));
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (200*dt*perme*ConstModulus/he^2)^(RKMethod) ));
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (200*dt*perme*ConstModulus/he^2/0.25)^RKMethod ));
            
            if ( implicit)
                AlphaStab = -0.65*perme*dt/he^2;

                if (AlphaStab < 0)
                    AlphaStab = 0;
                end
            end
        elseif ( all(ElementType == 'T6T6'))
            he = sqrt( sum( [GPInfo(el,:).Weight]) );
%             AlphaStab = 80*perme*dt/he^2;
            AlphaStab = 80*dt*perme/he^2*(1-exp(- (8E7*dt*perme/he^2)^(RKMethod) ));
            AlphaStab = 80*dt*perme/he^2*(1-exp(- (8E8*dt*perme/he^2)^(RKMethod) ));
            AlphaStab = 80*dt*perme/he^2*(1-exp(- (2000*dt*perme*ConstModulus/he^2/0.25)^(RKMethod) ));
        elseif ( all(ElementType == 'T6T3'))
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (6*dt*perme*ConstModulus/he^2/0.25)^(RKMethod) ));
        else
            disp(ElementType)
            error('this element does not exist. yet')
        end
        if (nargin == 11)
            AlphaStab = a/ConstModulus;%+b*perme*dt/he^2;
        end

        AlphaStab = -AlphaStab*AlphaStabM;


        Ms = GPInfo(el,ngp).Ms * AlphaStab;
        
        
        if ( all(ElementType=='T6T3') )
            Q2 = [Q'; zeros(3,12)];
            Ms = [Ms, zeros(3,3);
                -1,-1,0,2,0,0;
                 0, -1, -1, 0 ,2, 0;
                 -1, 0, -1, 0, 0, 2];
            Ce = [kke, Q, zeros(12,3); Q2, Ms];
            H = [H, zeros(3,3);
                zeros(3,6)];
            Ke = [zeros(12,18); 0*Q2, H];
        else
        
            Ce = [kke, Q; Q', Ms];
            Ke = [0*kke, 0*Q; 0*Q', H];
        end

        aux = GPInfo(el,ngp).IndexReorder;

        Ke = Ke(aux,aux);
        Ce = Ce(aux,aux);

        index = GPInfo(el).dofs;
        K(index,index) =  K(index,index) + Ke*GPInfo(el).Weight;
        C(index,index) =  C(index,index) + Ce*GPInfo(el).Weight;
    end
    
end

