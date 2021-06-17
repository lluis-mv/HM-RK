

syms xR_R real;
syms xR_V real;
syms xV_R real;
syms xV_V real;
syms xR real;
syms xV real;
syms R;

syms uR_r real;
syms uR_v real;
syms uV_r real;
syms uV_v real;

F = [xR_R, xR_V, 0;
    xV_R, xV_V, 0;
    0, 0, xR/R];

inv_F = inv(F);

J = F-eye(3)
j = eye(3)-inv_F

e = 0.5*(eye(3)-inv_F'*inv_F);


%% NOW REAL
syms R real;
syms Z real;
syms t real;
r = R*Z^2*exp(t)+R;
z = R*Z+2*t*Z*R*Z+Z;

F = [ diff(r,R), diff(r,Z), 0;
    diff(z,R), diff(z,Z), 0;
    0,0,r/R];

inv_F = inv(F)
J = F-eye(3)
j = eye(3)-inv_F

e = 0.5*(eye(3)-inv_F'*inv_F);

e_LMV = j(2,2);
for i = 1:3
    e_LMV = e_LMV - 0.5*( j(i,2)^2);
end