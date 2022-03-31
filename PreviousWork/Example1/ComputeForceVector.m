function [f, uDir, AllZero] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP)

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


f = zeros(3*nNodes, 1);
uDir = f;
AllZero = true;