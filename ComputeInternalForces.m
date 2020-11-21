
function [f] = ComputeInternalForces(Elements, GPInfo, X)


nElements = size(Elements, 1);

f = zeros(size(X));
m = [1;1;0];
for el = 1:nElements
    
    for gp = 1:size(GPInfo, 2)
               
        wP = GPInfo(el, gp).N * X( GPInfo(el,gp).dofsWP );
        
        index = GPInfo(el,gp).dofsU;
        
        f(index) = f(index) + GPInfo(el,gp).B'*( GPInfo(el,gp).StressNew([1,2,4]) - wP*m )*GPInfo(el,gp).Weight;

    end
end


