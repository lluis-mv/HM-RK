% Compute stress, stiffness....

function GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, U, C, consistent, RKMethod)

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

            GPInfo(el,gp) = EvaluateLaw( CP, GPInfo(el,gp), consistent, RKMethod);
        end

    end

else

    Idev = eye(3,3)-1/2*[ones(2,2), zeros(2,1); zeros(1,3)];
    II = 1/2*[ones(2,1); zeros(1,1)];
    for el = 1:nElem
        for gp = 1:size(GPInfo,2)

            nSystem = GPInfo(el,gp).dofsU;
            
            Uel = U(nSystem);
            Vol = U(GPInfo(el,gp).dofsVol);

            GPInfo(el,gp).StrainNew([1,2,4]) = GPInfo(el,gp).B*Uel;
            VolStrain = GPInfo(1,1).N*Vol;

            GPInfo(el,gp).AssumedStrainNew([1,2,4]) = Idev*GPInfo(el,gp).StrainNew([1,2,4]) + II*VolStrain;

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