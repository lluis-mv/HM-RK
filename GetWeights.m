function [al, be, w] = GetWeights(hydromechanical, element)

if ( nargin == 1)
    [al, be, w] = GetTheseWeights(hydromechanical);
else
    if ( hydromechanical)
        if ( all( element == 'T3T3') )
            [al, be, w] = GetTheseWeights(1);
        elseif (element(1)=='T')
            [al, be, w] = GetTheseWeights(3);
        elseif (all( element == 'Q8Q8'))
            [al, be, w] = GetTheseWeights(4);
        end
    else
        if ( all( element == 'T3T3') )
            [al, be, w] = GetTheseWeights(3);
        else
            [al, be, w] = GetTheseWeights(6);
        end
    end
end

function [al, be, w] = GetTheseWeights(nn)

if (nn == 1)
    al = 1/3;
    be = 1/3;
    w = 1;
elseif ( nn == 3)
    al = [1/6,1/6,2/3];
    be = [2/3,1/6,1/6];
    w = [1/3,1/3,1/3];
elseif ( nn == 6)

    wa = 0.054975871827661;
    wb = 0.1116907948390055;
    Na1 = 0.816847572980459;
    Nb1 = 0.108103018168070;
    Na2 = 0.091576213509771;
    Nb2 = 0.445948490915965;
    
    auxK = [Na2, Na2, wa;
        Na1, Na2, wa;
        Na2, Na1, wa;
        Nb2, Nb2, wb ;
        Nb1, Nb2, wb ;
        Nb2, Nb1, wb ];
    
    al = auxK(:,1)';
    be = auxK(:,2)';
    w = auxK(:,3)'/sum(auxK(:,3)');
elseif ( nn == 4)
    al = 1/sqrt(3)*[-1, 1, 1, -1];
    be = 1/sqrt(3)*[-1, -1, 1, 1];
    w = 1/4*[1,1,1,1];
elseif ( nn == 9)
    tt = sqrt(0.6);
    al = [-tt, -tt, -tt, 0, 0, 0, tt, tt, tt];
    be = [-tt, 0, tt, -tt, 0, tt, -tt, 0, tt];
    w = 1/81*[25, 40, 25, 40, 64, 40, 25, 40, 25]/4;
end
