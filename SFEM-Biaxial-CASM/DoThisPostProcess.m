function ThisInfo = DoThisPostProcess( t, Nodes, Elements, GPInfo, X, CP, PreviousInfo, finter)

if ( nargin < 8)
    finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
end


index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([index]-1)+2;

ThisInfo.t = t;
ThisInfo.F(1) = -sum(finter(dofs));
index = find( (Nodes(:,1) <= 1E-8) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([index]-1)+3;

ThisInfo.F(2) = X(dofs);

if (nargin >= 7)
    ThisInfo = [PreviousInfo, ThisInfo];
end
