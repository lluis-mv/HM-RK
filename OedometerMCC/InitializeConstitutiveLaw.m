
function [GPInfo] = InitializeConstitutiveLaw(GPInfo)

nElements = size(GPInfo,1);

for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el, gp).MCC = true;

        GPInfo(el, gp).StressNew =[10;10;10;0;0;0];
        GPInfo(el, gp).StressPrev = [10;10;10;0;0;0];

        GPInfo(el, gp).StrainNew = zeros(6,1);
        GPInfo(el, gp).StrainPrev = zeros(6,1);

        GPInfo(el, gp).HistoryNew = [500];
        GPInfo(el, gp).HistoryPrev = [500];
    end
end