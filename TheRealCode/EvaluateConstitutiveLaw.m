% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(GPInfo, U, C)

nElem = size(C,1);

for el = 1:nElem
    
    Ce = C(el,:);
    
    nSystem = [3*(Ce(1)-1)+[1,2],3*(Ce(2)-1)+[1,2], 3*(Ce(3)-1)+[1,2]];
    Uel = U(nSystem);
    
    GPInfo(el).StrainNew([1,2,4]) = GPInfo(el).B*Uel;
    
    GPInfo(el) = EvaluateLaw( GPInfo(el));
    
end



function GP = EvaluateLaw( GP)

if (GP.MCC)
    
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    
    X = [GP.StressPrev; GP.HistoryPrev];
    
    [Xnew, D] = ImplicitCamClay(X, DeltaStrain);
    
    GP.StressNew  = Xnew(1:6);
    GP.HistoryNew = Xnew(7);
    GP.D6 = D;
    GP.D = D([1,2,4], [1,2,4]);
    
else

    GP.StressNew = GP.StressPrev + GP.D6*(GP.StrainNew - GP.StrainPrev);
end