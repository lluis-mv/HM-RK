function [f, uDir, AllZero] = ComputeForceVector(t, Nodes, Elements, GPInfo, C, dofsPerNode)

if ( nargin == 5)
    dofsPerNode = 3;
end
if ( dofsPerNode ~= 3)
    error('please modify the code')
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(3*nNodes, 1);
uDir = f;
AllZero = true;
