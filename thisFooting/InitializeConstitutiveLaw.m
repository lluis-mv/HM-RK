function [GPInfo] = InitializeConstitutiveLaw(GPInfo, CP)

nElements = size(GPInfo,1);

addpath('../ModifiedCamClay');
[kappa, lambda, M, nu] = GetConstitutiveParameters();

Elastic = false;
if ( nargin == 2)
    Elastic = CP.Elastic;
end


for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el, gp).MCC = true;

        GPInfo(el, gp).StressNew =-[10;10;10;0;0;0];
        GPInfo(el, gp).StressPrev = -[10;10;10;0;0;0];

        GPInfo(el, gp).StrainNew = zeros(6,1);
        GPInfo(el, gp).StrainPrev = zeros(6,1);

        GPInfo(el, gp).HistoryNew = -[7.5];
        GPInfo(el, gp).HistoryPrev = -[7.5];
        
        p = -mean(GPInfo(el, gp).StressNew(1:3));
        if ( Elastic)
            K = p/kappa;
        else
            K = p/lambda;
        end
        
        GPInfo(el, gp).ConstrainedModulus =  3*K*(1-nu)/(1+nu);
    end
end