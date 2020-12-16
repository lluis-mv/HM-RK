
function [f] = ComputeInternalForces(Elements, GPInfo, X, HydroMechanical)

if (nargin == 3)
    error('Here you need the sign. ComputeInternal forces');
end

nElements = size(Elements, 1);

if ( HydroMechanical == true)
    sign = -1;
    Idev = eye(6);
else 
    sign = 1;
    Idev = eye(6)-1/3*[1,1,1,0,0,0]'*[1,1,1,0,0,0];
end

  


f = zeros(size(X));
m = [1;1;0];
for el = 1:nElements
    
    for gp = 1:size(GPInfo, 2)
               
        wP = GPInfo(el, gp).N * X( GPInfo(el,gp).dofsWP );
        
        index = GPInfo(el,gp).dofsU;
        ThisStress = Idev*GPInfo(el,gp).StressNew;
        f(index) = f(index) + GPInfo(el,gp).B'*( ThisStress([1,2,4]) + sign*wP*m )*GPInfo(el,gp).Weight;

    end
end



if ( HydroMechanical)
    return;
end

m6 = [1,1,1,0,0,0];

for el = 1:nElements
    
    for gp = 1:size(GPInfo, 2)
               
        wP = GPInfo(el, gp).N * X( GPInfo(el,gp).dofsWP );
        
        index = GPInfo(el,gp).dofsWP;
        
        f(index) = f(index) + GPInfo(el,gp).N' *m6*GPInfo(el,gp).StressNew*GPInfo(el,gp).Weight;
        f(index) = f(index) - 3*GPInfo(el,gp).N'* wP * GPInfo(el,gp).Weight;

    end
end


% i després em falta lo dels 6-3
theseDofs = [];
if ( length([GPInfo(1,1).dofsWP]) ~= length([GPInfo(1,1).dofsWPreal]) )
    
    for el = 1:nElements
        dofsWP = GPInfo(el,1).dofsWP;
        dofsReal = GPInfo(el,1).dofsWPreal(4:end);
        
        KK = 1/2*[1,1,0;
            0, 1, 1;
            1, 0,1];
        
        f(dofsReal) = -2*(KK*X(dofsWP)-X(dofsReal))*sum([GPInfo(el,:).Weight]);
        theseDofs = [theseDofs, dofsReal];
  
    end
    theseDofs = unique(sort(theseDofs));
    thisNorm = norm( f(theseDofs))
end

