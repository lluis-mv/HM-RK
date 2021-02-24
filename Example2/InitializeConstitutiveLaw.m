function [GPInfo] = InitializeConstitutiveLaw(GPInfo)

nElements = size(GPInfo,1);

addpath('../ModifiedCamClay');
[kappa, lambda, M, nu] = GetConstitutiveParameters();


for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el, gp).MCC = true;

        GPInfo(el, gp).StressNew =-[10;10;10;0;0;0];
        GPInfo(el, gp).StressPrev = -[10;10;10;0;0;0];

        GPInfo(el, gp).StrainNew = zeros(6,1);
        GPInfo(el, gp).StrainPrev = zeros(6,1);

        GPInfo(el, gp).HistoryNew = -[500];
        GPInfo(el, gp).HistoryPrev = -[500];
        
        p = -mean(GPInfo(el, gp).StressNew(1:3));
        K = p/kappa;
        
        GPInfo(el, gp).ConstrainedModulus =  3*K*(1-nu)/(1+nu);
    end
end





