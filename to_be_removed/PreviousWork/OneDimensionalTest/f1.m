function fTerm1 = f1(x1,x2,x3)
%F1
%    FTERM1 = F1(X1,X2,X3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-Mar-2022 15:55:08

t2 = x1.^2;
t3 = x1.^3;
t4 = x3.^2;
t5 = x3.^3;
t6 = x1.*(1.1e+1./3.0);
t7 = x3.*(2.2e+1./3.0);
t8 = t2.*(2.2e+1./3.0);
t10 = t3.*(3.3e+1./1.0e+1);
t9 = -t8;
fTerm1 = [t4.*(2.2e+1./3.0)-t5.*(3.3e+1./1.0e+1)+t6+t9+t10-x3.*(1.1e+1./3.0)+t2.*x3.*(1.1e+1./1.0e+1)-t4.*x1.*(1.1e+1./1.0e+1),t2.*(4.4e+1./3.0)-t3.*(3.3e+1./5.0)-t4.*(4.4e+1./3.0)+t5.*(3.3e+1./5.0)+t7-x1.*(2.2e+1./3.0)-t2.*x3.*(1.1e+1./5.0)+t4.*x1.*(1.1e+1./5.0),t4.*-2.2e+1+t5.*(6.6e+1./5.0)+t6+t7+t9+t10-x2.*1.1e+1-t2.*x2.*(1.1e+1./2.0)+t2.*x3.*(3.3e+1./5.0)+t4.*x1.*(9.9e+1./1.0e+1)-t4.*x2.*(3.3e+1./2.0)+x1.*x2.*(4.4e+1./3.0)-x1.*x3.*(4.4e+1./3.0)+x2.*x3.*(8.8e+1./3.0)-x1.*x2.*x3.*1.1e+1];
