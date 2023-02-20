function GP = FinalizeConstitutiveLaw(CP, GP)


for i = 1:size(GP,1)
    for gp = 1:size(GP,2)
        GP(i,gp).StrainPrev  = GP(i,gp).StrainNew;
        GP(i,gp).StressPrev  = GP(i,gp).StressNew;
        GP(i,gp).HistoryPrev = GP(i,gp).HistoryNew;
    end
end


if ( isfield( GP(1,1), 'AssumedStrainNew') )
    for i = 1:size(GP,1)
        for gp = 1:size(GP,2)
            GP(i,gp).AssumedStrainPrev  = GP(i,gp).AssumedStrainNew;
        end
    end
end





if ( GP(1,1).MCC )
    [kappa, lambda, M, nu, n, r, m] = GetConstitutiveParametersCASM();
    
    for i = 1:size(GP,1)
        for gp = 1:size(GP,2)
            p = -mean(GP(i, gp).StressNew(1:3));
            
            if ( CP.Elastic)
                K = p/kappa;
            else
                K = p/lambda;
            end
            
            GP(i, gp).ConstrainedModulus =  3*K*(1-nu)/(1+nu);
        end
    end
end