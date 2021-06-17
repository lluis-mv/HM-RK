function s = solution(XX, TT)
s = XX.^6+(XX).*(XX-1);
TT
if ( TT > 1)
    return;
end
%function pw = EvalConsolidation(XX, TT)
pw = 0*XX;
nTermsZero = 0;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:5000
            aux = pi/2*(2*m+1);
            term = 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT);
            if ( abs(term) < 1e-12)
                nTermsZero = nTermsZero + 1;
            else
                nTermsZero = 0;
            end
            pw(i,j) = pw(i,j)+term ;
            if ( nTermsZero > 20)
                break;
            end
        end
    end
end

s = pw;