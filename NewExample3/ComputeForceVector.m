function [f, uDir, AllZero, F] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP, dofsPerNode)

if ( nargin == 5)
    dofsPerNode = 3;
end

if (dofsPerNode~=3)
    error('This has not been yet implemented')
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(3*nNodes, 1);



index = find( (Nodes(:,1) <= 1) &  abs(Nodes(:,2)) < 1E-8);
uDir = 0*f;
uDir( 3*(index-1)+2) = -0.05;

AllZero = false;

F = 0*f;