function ThisInfo = DoThisPostProcess( t, Nodes, Elements, GPInfo, X, CP, PreviousInfo, finter, dofsPerNode)

if ( nargin < 8)
finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
end

if (  length(finter) == 0)
    finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
end
if ( nargin < 9)
    dofsPerNode = 3;
end

if ( dofsPerNode ~= 3)
    error('You need to implement this')
end

index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([index]-1)+2;

ThisInfo.t = t;
ThisInfo.F(1) = -sum(finter(dofs));
index = find( (Nodes(:,1) <= 1E-8) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([index]-1)+3;

ThisInfo.F(2) = X(dofs);
ThisInfo.F(3) = -X(dofs-1); 

if (nargin == 7)
    ThisInfo = [PreviousInfo, ThisInfo];
end
