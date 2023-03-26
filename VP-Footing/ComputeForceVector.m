function [f, uDir, AllZero, F] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP, dofsPerNode)

if (nargin == 5)
    dofsPerNode = 3;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(dofsPerNode*nNodes, 1);



index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
uDir = 0*f;
uDir( dofsPerNode*(index-1)+2) = -0.10/3600.0;

AllZero = false;

F = 0*f;
