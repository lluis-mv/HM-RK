function [] = main2()





function k = ComputeKProblem(k,y,a,b,dt)

f = @(x) ComputeKError(x,y,a,b,dt);
options = optimoptions('fsolve', 'FunctionTolerance', 1E-14, 'StepTolerance', 1E-14);
k = fsolve(f, k,options);
k = fsolve(f, k,options);
k = fsolve(f, k,options);

function err = ComputeKError(k,y,a,b,dt)


kN = zeros( 3, length(b));


    for i = 1:length(b)
        yS = y;
        for j = 1:length(b)
            yS = yS + dt*a(i,j)*k(1,j);
        end
        kN(:,i) = sourceTerm(yS);
    end

err = k-kN;
% err = norm(k-kN);



function dy = sourceTerm(y)
dy = y;



function [a, b] = GetImplRungeKutta(RK)

if ( RK == -1)
    a = 1;
    b = 1;
elseif ( RK == -2)
    a = [0,0;
        1/2, 1/2];
    b = [1/2, 1/2];
elseif (RK == -3)
    a = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
        5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
        5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
    b = [5/18, 4/9, 5/18];
end