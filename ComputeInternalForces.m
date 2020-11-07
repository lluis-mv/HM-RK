
function [f] = ComputeInternalForces(Elements, GPInfo, X)


nElements = size(Elements, 1);

f = zeros(size(X));
m = [1;1;0];
for el = 1:nElements
    
    for gp = 1:size(GPInfo, 2)
%         ind = Elements(el,:);
%         index = [];
%         for ii = 1:length(ind)
%             index = [ index, (ind(ii)-1)*3 + [1,2] ];
%         end

%         PW = X( 3*(ind-1)+3);

%         pw = 1/3*(PW(1)+PW(2)+PW(3));

        %f(index) = f(index) + GPInfo(el).B'*( GPInfo(el).StressNew([1,2,4]) + m*pw)*GPInfo(el).Weight;
        
        index = GPInfo(el,gp).dofsU;
        
        f(index) = f(index) + GPInfo(el,gp).B'*( GPInfo(el,gp).StressNew([1,2,4]) )*GPInfo(el,gp).Weight;
        

    end
end


disp('he de fer aquesta funció, saps....')