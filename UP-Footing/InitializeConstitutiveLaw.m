function [GPInfo] = InitializeConstitutiveLaw(GPInfo)

nElements = size(GPInfo,1);

for el = 1:nElements
    for gp = 1:size(GPInfo,2)
        GPInfo(el, gp).VonMises = true;
    end
end