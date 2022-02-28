% Ensamble elemtal matrices to create C and K
% Compute the stabilization factor

function [C, K] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM)

if (nargin == 8)
    AlphaStabM = 1;
end

if ( CP.HydroMechanical)
    [C, K] = EnsambleHydroMechanicalProblem(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM);
else
    [C, K] = EnsambleUPProblem(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM);
end


function [C, K] = EnsambleHydroMechanicalProblem(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM)


nNodes = size(Nodes, 1);
nElements = size(Elements, 1);
nSystem = 3*nNodes;

C = sparse(nSystem, nSystem);
K = C;


one = [1,1,0]';

perme = CP.k;




for el = 1:nElements
    
    for ngp = 1:size(GPInfo,2)
        
        
        ConstModulus=  GPInfo(el,ngp).ConstrainedModulus;
        
        kke = GPInfo(el,ngp).B'*GPInfo(el,ngp).D*GPInfo(el,ngp).B;
        Q = -GPInfo(el,ngp).B'*one * GPInfo(el,ngp).N;
        H = GPInfo(el,ngp).dN_dX'*perme*GPInfo(el,ngp).dN_dX;
    
        
        he = sqrt( sum([GPInfo(el,:).Weight]));
   
        if (all(ElementType == 'T3T3'))
            term = exp(- (200*dt*perme*ConstModulus/he^2/0.25)^RKMethod);
            AlphaStab = 8*dt*perme/he^2*(1-term) + ConstModulus/1000000*term;
        elseif ( all(ElementType == 'T6T6'))
            term = exp(- (2000*dt*perme*ConstModulus/he^2/0.25)^(RKMethod) );
            AlphaStab = +80*dt*perme/he^2*(1-term) +ConstModulus/1000000*term;
        elseif ( all(ElementType == 'T6T3'))
            AlphaStab = 8*dt*perme/he^2*(1-exp(- (6*dt*perme*ConstModulus/he^2/0.25)^(RKMethod) ));
        else
            AlphaStab = 0;
        end
        
       

        if ( implicit)
            if ( all(ElementType == 'T3T3'))
                AlphaStab = 2/ConstModulus - 12*dt*perme/he^2;
                AlphaStab = max(0, AlphaStab);
            else
                AlphaStab = 0;
            end
        end
        

        
        
        if ( length(AlphaStabM) == 1)
            AlphaStab = -AlphaStab*AlphaStabM;
        end
        
        
        
        if ( length(AlphaStabM) == 1 && AlphaStabM(1) < 0)
            if (all(ElementType == 'T3T3'))
                term = exp(- (200*dt*perme*ConstModulus/he^2/0.25)^RKMethod);
                AlphaStab = abs(AlphaStabM)*8*dt*perme/he^2 + ConstModulus/1000000*term;
            elseif ( all(ElementType == 'T6T6'))
                term = exp(- (2000*dt*perme*ConstModulus/he^2/0.25)^(RKMethod) );
                AlphaStab = abs(AlphaStabM)*80*dt*perme/he^2 +ConstModulus/1000000*term;
            elseif ( all(ElementType == 'T6T3'))
                AlphaStab = abs(AlphaStabM)*8*dt*perme/he^2;
            end
            AlphaStab = -AlphaStab;
        end
     
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
        elseif ( all(ElementType =='Q8Q4') )
            Q2 = [Q'; zeros(4,16)];
            Ms = [0.0*Ms, zeros(4,4);
                 -1,-1, 0, 0,  2,0,0,0;
                  0,-1,-1, 0,  0,2,0,0;
                  0, 0,-1,-1,  0,0,2,0;
                 -1, 0, 0,-1,  0,0,0,2];
            Ce = [kke, Q, zeros(16,4); Q2, Ms];
            H = [H, zeros(4,4);
                zeros(4,8)];
            Ke = [zeros(16,24); 0*Q2, H];
        else
        
            Ce = [kke, Q; Q', Ms];
            Ke = [0*kke, 0*Q; 0*Q', H];
        end

        aux = GPInfo(el,ngp).IndexReorder;

        Ke = Ke(aux,aux);
        Ce = Ce(aux,aux);

        index = GPInfo(el,ngp).dofs;
        K(index,index) =  K(index,index) + Ke*GPInfo(el,ngp).Weight;
        C(index,index) =  C(index,index) + Ce*GPInfo(el,ngp).Weight;
    end
end


function [C, K] = EnsambleUPProblem(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, implicit, AlphaStabM)


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
        

        AlphaStab = -800*AlphaStabM;
        AlphaStab = 0;

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

K = 0*C;
