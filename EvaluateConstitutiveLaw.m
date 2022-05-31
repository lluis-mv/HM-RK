% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, U, C, consistent, RKMethod)

if ( nargin == 5)
    RKMethod = 0;
end
nElem = size(C,1);

for el = 1:nElem
    for gp = 1:size(GPInfo,2)
    
        nSystem = GPInfo(el,gp).dofsU;
        Uel = U(nSystem);
    
        GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;
    
        GPInfo(el,gp) = EvaluateLaw( CP, GPInfo(el,gp), consistent, RKMethod);
    end
    
end

