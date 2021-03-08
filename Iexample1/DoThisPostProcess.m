function ThisInfo = DoThisPostProcess( t, Nodes, Elements, GPInfo, X, CP, PreviousInfo)

ThisInfo = [];
if (nargin == 7)
    ThisInfo = [PreviousInfo, ThisInfo];
end
