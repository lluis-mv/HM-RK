function [GPInfo] = InitializeConstitutiveLaw(CP, GPInfo)

nElements = size(GPInfo,1);

addpath('../ModifiedCamClay');
[kappa, lambda, M, nu] = GetConstitutiveParameters();


Elastic = CP.Elastic;

number = 1;
if ( length([GPInfo(1,1).dofsWP]) == length([GPInfo(1,1).dofsWPreal]) )
    if ( length([GPInfo(1,1).dofsWPreal]) == 3 )
        number = 20;
    else
        number = 120;
    end
end

for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el, gp).MCC = CP.MCC;

        GPInfo(el, gp).StressNew =-[10;10;10;0;0;0];
        GPInfo(el, gp).StressPrev = -[10;10;10;0;0;0];

        GPInfo(el, gp).StrainNew = zeros(6,1);
        GPInfo(el, gp).StrainPrev = zeros(6,1);

        GPInfo(el, gp).HistoryNew = -11.0;
        GPInfo(el, gp).HistoryPrev = -11.0;
        
        p = -mean(GPInfo(el, gp).StressNew(1:3));
        if ( Elastic)
            K = p/kappa;
        else
            K = p/lambda*number;
        end
        
        GPInfo(el, gp).ConstrainedModulus =  3*K*(1-nu)/(1+nu);
    end
end