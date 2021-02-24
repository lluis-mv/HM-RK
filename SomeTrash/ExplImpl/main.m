function [] = main()



RKMethods = [1:8];
RKRef = 8;
RKMethods = [-1,-2,-3,8];
RKMethods = [-1,-2,-3,-4,-5,-6,8];
NSTEPS = 2.^[0:8];

nS = 1;
for nSteps = NSTEPS
    
    rk = 1;
    for RK = RKMethods
        
        y = IntegrateThis(RK, nSteps);
        YY(rk, nS) = y;
        
        if ( RK == RKRef)
            YYadim = y;
        end
        
        rk = rk+1;
    end
    figure(1)
    hold off
    for i = 1:length(RKMethods)
        loglog(NSTEPS(1:nS), YY(i,1:nS), '*-.')
        hold on
    end
    
    figure(2)
    hold off
    for i = 1:length(RKMethods)
        loglog(NSTEPS(1:nS), abs(YY(i,1:nS)-YYadim), '*-.', 'DisplayName', num2str(RKMethods(i)) )
        hold on
    end
    legend('location', 'best')
    
    nS = nS+1;
end

function y = IntegrateThis(RK, nSteps)

if ( RK > 0)
    y = IntegrateExpl(RK, nSteps);
    return;
end

dt = 1/nSteps;
[a, b] = GetImplRungeKutta(RK);
k = rand(  1, length(b));
y = 1;
for step = 1:nSteps
    k = ComputeKProblem(k, y, a, b, dt);
    for i = 1:length(b)
        y = y+dt*b(i)*k(1,i);
    end
end



function k = ComputeKProblem(k,y,a,b,dt)

f = @(x) ComputeKError(x,y,a,b,dt);
options = optimoptions('fsolve', 'FunctionTolerance', 1E-14, 'StepTolerance', 1E-14);
k = fsolve(f, k,options);
k = fsolve(f, k,options);
k = fsolve(f, k,options);

function err = ComputeKError(k,y,a,b,dt)


kN = zeros(  1, length(b));


    for i = 1:length(b)
        yS = y;
        for j = 1:length(b)
            yS = yS + dt*a(i,j)*k(1,j);
        end
        kN(1,i) = sourceTerm(yS);
    end

err = k-kN;
% err = norm(k-kN);

function y = IntegrateExpl(RK, nSteps)

[a, b] = GetRungeKutta(RK);
k = zeros(  1, length(b));
y = 1;
dt = 1/nSteps;
for step = 1:nSteps
    for i = 1:length(b)
        yS = y;
        for j = 1:i-1
            yS = yS + dt*a(i,j)*k(1,j);
        end
        k(1,i) = sourceTerm(yS);
    end
    
    for i = 1:length(b)
        y = y+dt*b(i)*k(1,i);
    end
    
end


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
    b = [-5/6, 8/3, -5/6];
elseif (RK == -4)
    a = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
        5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
        5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
    b = [5/18, 4/9, 5/18];
elseif (RK == -5)
    a = [1/6, -1/6, 0;
        1/6, 1/3, 0;
        1/6, 5/6, 0];
    b = [1/6, 2/3, 1/6];
elseif (RK == -6)
    a = [1/9, (-1-sqrt(6))/18, (-1+sqrt(6))/18; 
        1/9, 11/45+(7*sqrt(6))/360, 11/45-(43*sqrt(6))/360; 
        1/9, 11/45+(43*sqrt(6))/360, 11/45-(7*sqrt(6))/360]; 
%         b = [1/9, 4/9+sqrt(6)/36, 4/5-sqrt(6)/36];
        b = [1/9, 4/9+sqrt(6)/36, 4/9-sqrt(6)/36];
end