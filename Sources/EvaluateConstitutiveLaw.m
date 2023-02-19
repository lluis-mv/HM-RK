% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, U, C, consistent, RKMethod, DT)

if ( nargin == 5)
    RKMethod = 0;
end

nElem = size(C,1);

if ( isfield(GPInfo(1,1), 'AssumedStrainNew') == false)

    for el = 1:nElem
        for gp = 1:size(GPInfo,2)

            nSystem = GPInfo(el,gp).dofsU;
            Uel = U(nSystem);

            GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;

            GPInfo(el,gp) = EvaluateLaw( CP, GPInfo(el,gp), consistent, RKMethod, DT);
        end

    end

else

    Idev = eye(6,6)-1/3*[ones(3,3), zeros(3,3); zeros(3,6)];
    II = 1/3*[ones(3,1); zeros(3,1)];
    for el = 1:nElem
        for gp = 1:size(GPInfo,2)

            nSystem = GPInfo(el,gp).dofsU;
            
            Uel = U(nSystem);
            Vol = U(GPInfo(el,gp).dofsVol);

            GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;
            VolStrain = GPInfo(1,1).N*Vol;

            GPInfo(el,gp).AssumedStrainNew = Idev*GPInfo(el,gp).StrainNew + II*VolStrain;

            Aux1 = GPInfo(el,gp).StrainNew;
            Aux2 = GPInfo(el,gp).StrainPrev;

            GPInfo(el,gp).StrainNew = GPInfo(el,gp).AssumedStrainNew;
            GPInfo(el,gp).StrainPrev = GPInfo(el,gp).AssumedStrainPrev;
            GPInfo(el,gp) = EvaluateLaw( CP, GPInfo(el,gp), consistent, RKMethod);


            GPInfo(el,gp).StrainNew = Aux1;
            GPInfo(el,gp).StrainPrev = Aux2;
        end
    end


end