
function [GPInfo] = InitializeConstitutiveLaw(GPInfo)

nElements = length(GPInfo);

for el = 1:nElements
    GPInfo(el).MCC = true;
    
    GPInfo(el).StressNew =[10;10;10;0;0;0];
    GPInfo(el).StressPrev = [10;10;10;0;0;0];
    
    GPInfo(el).StrainNew = zeros(6,1);
    GPInfo(el).StrainPrev = zeros(6,1);
    
    GPInfo(el).HistoryNew = [500];
    GPInfo(el).HistoryPrev = [500];
end