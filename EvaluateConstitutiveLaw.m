% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(GPInfo, U, C, consistent)

nElem = size(C,1);

for el = 1:nElem
    for gp = 1:size(GPInfo,2)
    
        Ce = C(el,:);
    
        nSystem = GPInfo(el,gp).dofsU;
        Uel = U(nSystem);
    
        GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;
    
        GPInfo(el,gp) = EvaluateLaw( GPInfo(el,gp), consistent);
    end
    
end



function GP = EvaluateLaw( GP, implicit)

if (GP.MCC)
    
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    
    X = [GP.StressPrev; GP.HistoryPrev];
    
    %[Xnew, D] = ImplicitCamClay(X, DeltaStrain);
    [Xnew, Dconsist, D] = ExplicitCamClay(X, DeltaStrain, -1);
    if (implicit)
        D = Dconsist;
    end
    
    GP.StressNew  = Xnew(1:6);
    GP.HistoryNew = Xnew(7);
    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
    
else

    GP.StressNew = GP.StressPrev + GP.D6*(GP.StrainNew - GP.StrainPrev);
end