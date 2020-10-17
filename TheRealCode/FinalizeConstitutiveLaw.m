function GP = FinalizeConstitutiveLaw(GP)


for i = 1:length(GP)
    GP(i).StrainPrev  = GP(i).StrainNew;
    GP(i).StressPrev  = GP(i).StressNew;
    GP(i).HistoryPrev = GP(i).HistoryNew;
end
