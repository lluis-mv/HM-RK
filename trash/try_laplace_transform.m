function []=try_laplace_transform
close all;



gamma = 10000.00;
gamma = 0;
TT = linspace(1e-5,5e-4,100);
p = 0*TT;

for Inff = [1e4,1e5,1e6,1e7,1e8]
    smallNumber = 1e-9;
    myEps = 1e-11;
    for ii = 1:length(TT)
        fun = @(s) source(s, TT(ii));
        p(ii) = integral(fun, gamma-Inff*1i,gamma-smallNumber*1i)/2/pi/1i;
        p(ii) = p(ii)+ integral(fun, gamma-smallNumber*1i,gamma-myEps*1i)/2/pi/1i;
        p(ii) = integral(fun, gamma-myEps*1i, gamma+myEps*1i)/2/pi/1i;
        p(ii) = p(ii)+ integral(fun, gamma+myEps*1i,gamma+smallNumber*1i)/2/pi/1i;
        p(ii) = p(ii) + integral(fun, gamma+smallNumber*1i,gamma+Inff*1i)/2/pi/1i;
        
    end
    plot(TT,p);
    hold on
    legend('1e4','1e5','1e6','1e7','1e8','1e9', 'location','best')
    pause(0.0001)
end


function p = source(s, t)

% syms rho_w positive
% syms rho_s positive
% syms n positive
% syms Cf positive
% syms Cs positive
% syms m_v positive
% syms k positive
% syms s
% syms t
% syms x
x = 0.2;


n = 0.4;
alfa = 1;
rho_w = 1;
rho_s = 2.65;

k = 0.00005;

m_v = 2e-10;
Cf = 5e-10;
Cs = 0;

S_p = n * Cf + (alfa -n)*Cs;

S_p = 1/(5E6);
m_v = 1/5E6;

A = -rho_w * S_p * s.^2;
B = (rho_w-rho_s) * ( 1-n)*m_v * s.^2;
C = -(rho_w * s + n/k) .* (S_p * s/n);
D = ((rho_w * s + n/k)*(1-n)/n + n / k )*m_v.*s;

X = B + C - D;
Y = B.*C-A.*D;


alpha1_2 =  sqrt(X.^2-4*Y) / 2 - X/2;
alpha2_2 =  -sqrt(X.^2-4*Y) / 2 - X/2;
E1 = -(1./s).*(   (alpha2_2+C)./(alpha1_2-alpha2_2));
F1 = (1./s).*(   (alpha1_2+C)./(alpha1_2-alpha2_2));

alpha1 = sqrt(alpha1_2);
alpha2 = sqrt(alpha2_2);
p = E1.*exp(-alpha1.*x) + F1.*exp(-alpha2.*x);
p = exp(s*t).*p;
disp('hola')

