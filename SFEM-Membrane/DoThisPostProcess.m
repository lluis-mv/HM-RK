function ThisInfo = DoThisPostProcess( t, Nodes, Elements, GPInfo, X, CP, PreviousInfo, finter)

if ( nargin < 8)
    finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);
end


index = find( Nodes(:,1) == max(Nodes(:,1)) &  Nodes(:,2) == max(Nodes(:,2)));
dofs = 3*([index]-1)+2;

ThisInfo.t = t;
ThisInfo.F(1) = -sum(X(dofs));
dofs = 3*([index]-1)+3;
ThisInfo.F(2) = X(dofs);
if (nargin >= 7)
    ThisInfo = [PreviousInfo, ThisInfo];
end
