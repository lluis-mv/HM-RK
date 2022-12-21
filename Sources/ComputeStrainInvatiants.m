function GPInfo = ComputeStrainInvatiants(GPInfo)

for i = 1:size(GPInfo,1)
    for j = 1:size(GPInfo, 2)
        Strain = GPInfo(i,j).StrainNew;
        eVol = sum(Strain(1:3));
        eDev = 0;
        for k = 1:3
            eDev = eDev + (Strain(k)-eVol/3)^2;
        end
        for k = 4:6
            eDev = eDev + 2*(1/2*Strain(k))^2;
        end
        GPInfo(i,j).StrainVol = eVol;
        GPInfo(i,j).StrainDev = sqrt(eDev);
    end
end