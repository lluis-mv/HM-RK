% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(GPInfo, U, C, consistent, RKMethod)

if ( nargin == 4)
    RKMethod = 0;
end
nElem = size(C,1);

for el = 1:nElem
    for gp = 1:size(GPInfo,2)
    
        nSystem = GPInfo(el,gp).dofsU;
        Uel = U(nSystem);
    
        GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;
    
        GPInfo(el,gp) = EvaluateLaw( GPInfo(el,gp), consistent, RKMethod);
    end
    
end



function GP = EvaluateLaw( GP, consistent, RKMethod)

if (GP.MCC)
    
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    X = [GP.StressPrev; GP.HistoryPrev];
    
    if (consistent)
        [Xnew, D, ~] = ExplicitCamClay2(X, DeltaStrain, -1);
    else
        if ( RKMethod == 1)
            RKMethod = RKMethod+1;
        end
        [Xnew, ~, D] = ExplicitCamClay2(X, DeltaStrain, RKMethod, false);
    end
    
    
    GP.StressNew  = Xnew(1:6);
    GP.HistoryNew = Xnew(7);
    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
elseif (GP.VonMises)
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    X = [GP.StressPrev];
    
    [Xnew, Dconsist, D] = ExplicitVonMises(X, DeltaStrain, -1);
    if (consistent)
        D = Dconsist;
    end
    
    GP.StressNew  = Xnew(1:6);

    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
    
else
    GP.StressNew = GP.StressPrev + GP.D6*(GP.StrainNew - GP.StrainPrev);
end
