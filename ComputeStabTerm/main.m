syms alfa real
syms beta real


N =  [ 1 - alfa - beta; alfa;  beta];

% N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
%     alfa*(2*alfa-1);
%     beta*(2*beta-1);
%     4*alfa*(1-alfa-beta);
%     4*alfa*beta;
%     4*beta*(1-alfa-beta)];

Proj = int( N, beta, 0, 1-alfa);
Proj = int( Proj, alfa, 0, 1);

Area = int( 1, beta, 0, 1-alfa);
Area = int( Area, alfa, 0, 1);

Proj = Proj/Area;


Ms = int( (Proj-N)*(Proj-N)', beta, 0, 1-alfa)
Ms = int( Ms, alfa, 0, 1)



N = 1/4*[(1-alfa)*(1-beta); (1+alfa)*(1-beta); (1+alfa)*(1+beta); (1-alfa)*(1+beta)]

Proj = int( N, beta, -1, 1);
Proj = int( Proj, alfa, -1, 1);

Area = int( 1, beta, -1, 1);
Area = int( Area, alfa, -1, 1);

Proj = Proj/Area;


Ms = int( (Proj-N)*(Proj-N)', beta, -1, 1)
Ms = int( Ms, alfa, -1, 1)

% Now check this stupid integration
for nn = [4, 9]
    if ( nn == 4)
        al = 1/sqrt(3)*[-1, 1, 1, -1];
        be = 1/sqrt(3)*[-1, -1, 1, 1];
        w = 1/4*[1,1,1,1];
    elseif ( nn == 9)
        tt = sqrt(0.6);
        al = [-tt, -tt, -tt, 0, 0, 0, tt, tt, tt];
        be = [-tt, 0, tt, -tt, 0, tt, -tt, 0, tt];
        w = 1/81*[25, 40, 25, 40, 64, 40, 25, 40, 25];
    end
    Nint = 0*N;
    
    for i = 1:length(al)
        thisAlfa = al(i);
        thisBeta = be(i);
        weigth = w(i);
        
        term = subs( N, alfa, thisAlfa);
        term = subs(term, beta, thisBeta);
        Nint = Nint + term*weigth;
        Nint = simplify(Nint);
    end
end
        