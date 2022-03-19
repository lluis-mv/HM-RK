function [u, pw] = EvaluateSticklePastor( x, t, L, cv, M)

if ( nargin == 0)
    figure(1); clf;
    xx = linspace(0,1,10);
    for j = 1:length(xx)
        tt = linspace(0,4);
        pw = tt*0;
        u = tt*0;
        for i = 1:length(tt)
            [u(i), pw(i)] =EvaluateSticklePastor(xx(j), tt(i), 1, 1,1);
        end
        subplot(1,2,1)
        plot(tt, u)
        hold on
        subplot(1,2,2)
        plot(tt, pw)
        hold on
    end
    return;
end

tc = 1;
q =  1;

pw = 0;
f =  q*t/tc*(1 - heaviside(t-tc)) + q*heaviside(t-tc);
u = f*(L-x)/M;
for n = 0:200
    
    Nn = cv*( (1+2*n)*pi/2/L)^2;
    
    if ( t < tc)
        h = 1 - exp(-Nn*t);
    else
        h = exp(-Nn*(t-tc))-exp(-Nn*t);
    end
    
    pn = 4*q*h/pi/(2*n+1)/tc/Nn;
    pn2 = pn*sin( (1+2*n)*pi*x/2/L);
    pw = pw+pn2;
    
    
    
    u = u - 1 / M * 2*L/(1+2*n)/pi*pn * cos( (1+2*n)*pi*x/2/L);
end