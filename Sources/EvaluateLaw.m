

function GP = EvaluateLaw(CP, GP, consistent, RKMethod)

if (GP.MCC)
    
    DeltaStrain = GP.StrainNew-GP.StrainPrev;
    
    X = [GP.StressPrev; GP.HistoryPrev];
    
    if (consistent)
        if ( CP.Elastic)
            [Xnew, D, ~] = ExplicitCamClayE(X, DeltaStrain, -1);
        elseif ( GP.MCC == 1)
            [Xnew, D, ~] = Hashiguchi3(X, DeltaStrain, CP.RK, true);
        elseif ( GP.MCC == 2)
            [Xnew, D, ~] = ExplicitCamClay2(X, DeltaStrain, -1);
        elseif ( GP.MCC == 20)
            [Xnew, D] = ImplicitCASM(X, DeltaStrain);
        end
    else
        if ( RKMethod == 1)
            RKMethod = RKMethod+1;
        end
        if ( CP.Elastic)
            [Xnew, ~, D] = ExplicitCamClayE(X, DeltaStrain, RKMethod, false);
        else
            [Xnew, ~, D] = ExplicitCamClay2(X, DeltaStrain, RKMethod, false);
        end
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
