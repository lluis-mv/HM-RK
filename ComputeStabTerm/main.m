syms alfa real
syms beta real

syms E real
syms nu real

N =  [ 1 - alfa - beta; alfa;  beta];

for i = 1:3
    for j = 1:3
        if ( i==j)
            De(i,j) = (1-nu);
        else
            De(i,j) = nu;
        end
    end
    De(i+3, i+3) = (1-2*nu)/2;
end

De = E/(1+nu)/(1-2*nu) * De;


D = De([1,2,4], [1,2,4]);


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


Ms = int( (Proj-N)*(Proj-N)', beta, -1, 1);
Ms = int( Ms, alfa, -1, 1)

M = int(N*N', beta, -1, 1);
M = int(M, alfa, -1,1)


gradN = [diff(N, alfa), diff(N, beta)];
dN_dX = gradN';


H = int( gradN*gradN',  beta, -1, 1);
H = int( H, alfa, -1, 1)


    B = [];
        for i = 1:4
            b = [dN_dX(1,i), 0; 0, dN_dX(2,i); dN_dX(2,i), dN_dX(1,i)]; %% geotechnical engineering
            B = [B, b];
        end




return;

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
        