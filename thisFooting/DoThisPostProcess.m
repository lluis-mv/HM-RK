function ThisInfo = DoThisPostProcess( t, Nodes, Elements, GPInfo, X, CP, PreviousInfo)

finter = ComputeInternalForces( Elements, GPInfo, X, CP.HydroMechanical);


index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
dofs = 3*([index]-1)+2;

ThisInfo.t = t;
ThisInfo.F = -sum(finter(dofs));


if (nargin == 7)
    ThisInfo = [PreviousInfo, ThisInfo];
end


figure(232)
plot([ThisInfo.t], [ThisInfo.F], '*-.')
hold on
drawnow;