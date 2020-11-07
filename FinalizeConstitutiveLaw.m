function GP = FinalizeConstitutiveLaw(GP)


for i = 1:size(GP,1)
    for gp = 1:size(GP,2)
        GP(i,gp).StrainPrev  = GP(i,gp).StrainNew;
        GP(i,gp).StressPrev  = GP(i,gp).StressNew;
        GP(i,gp).HistoryPrev = GP(i,gp).HistoryNew;
    end
end
