function [u,v,out3,x,y] = solution(x,y)
%SOLUTION
%    [U,V,OUT3,X,Y] = SOLUTION(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    11-Mar-2022 09:34:43

u = 0.0;
if nargout > 1
    v = (y.^4.*(y-1.0).^2)./1.0e+2;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    x = x;
end
if nargout > 4
    y = y;
end
