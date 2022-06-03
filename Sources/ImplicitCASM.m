function [Xnew, Dc] =ImplicitCASM(X, DeltaStrain)
% I guess I do solid mechanics sign convention. Not sure

global tolerance
tolerance = 5E-13;

[Xnew, Dc, plastic] = ReturnMapping(X, DeltaStrain);
Dc = Dc(1:6,1:6);


% % return;
% % [D, De, H] = ComputeContinuousMatrix(Xnew(1:6), Xnew(7), plastic);
% % if ( plastic)
% % hola = 1;
% % end


function [XNew, D, plastic] = ReturnMapping( X, DeltaStrain)

global tolerance

plastic = false;

[kappa, lambda, M, nu, n, r, m] = GetConstitutiveParametersCASM();

mIdentity = reshape([1.0,1.0,1.0,0.0,0.0, 0.0], [6,1]);

ThisDev = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0],[6,6]);


DeltaStrainVol = sum(DeltaStrain(1:3));
DeltaStrainDev = DeltaStrain - DeltaStrainVol/3*mIdentity;

p0 = sum(X(1:3))/3;
s0 = X(1:6)-p0*mIdentity;
pc0 = X(7);




p = p0 * exp( -DeltaStrainVol/kappa);
pc = pc0;

K = (-p)/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);

s = s0 + 2*G*ThisDev * DeltaStrainDev;

J2sqrt = 0;
for i = 1:6
    number = 1;
    if i > 3
        number = 2;
    end
    J2sqrt = J2sqrt + number*s(i)^2;
end
J2sqrt = sqrt(0.5*J2sqrt);
yield = EvaluateYieldSurface(p, pc, J2sqrt, M, n, r);


if (  yield < tolerance)
    XNew = [p*mIdentity+s; pc];
    
    dPdeVol = -p0 * exp( -DeltaStrainVol/kappa)/kappa;
    D = zeros(6,6);
    D(1:3,1:3) = dPdeVol;
    D = D + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3); zeros(3,3), 0.5*eye(3,3)];
    
    dGdEvol = 3*(1-2*nu)/2/(1+nu)*p0*exp(-DeltaStrainVol/kappa)/kappa^2;
    
    D(:,1:3) = D(:,1:3)  + 2*ThisDev*DeltaStrainDev*dGdEvol;
    D = [D; zeros(1,6)];
    

    return;
    
end

plastic = true;

gamma = 0;
iter = 0;
Jacobian = zeros(9,9);


while (true)
 
    [residual] = ComputeResidual(p, pc, gamma, s,  p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
    
    normResidual = norm(residual);
    %disp(['iter ', num2str(iter), ' RESI ', num2str(normResidual)])
    
    if ( normResidual < tolerance && iter > 1)
        break
    end
    
    if ( iter > 50)
        if ( normResidual > 100*tolerance)
            disp('maxIter')
        end
        break;
    end
    
    Jacobian = EvaluateJacobian(p, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
    
    
    
    if ( rcond(Jacobian) < 1E-15) || ( isnan(rcond(Jacobian) ))
        hola = 1;
    end
    
    dX = -Jacobian\residual;
    dX = real(dX);
    if (norm(dX) > 1)
        dX = dX/norm(dX);
    end
    
    
    p = p + dX(1);
    pc = pc+dX(2);
    gamma = gamma+dX(3);
    s = s + dX(4:9);
    
    if ( any(isnan(dX)))
        break;
    end
    iter = iter+1;
    
end
K = (-p)/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);

XNew = [p*mIdentity+s; pc];


% Now we have to obtain the stiffness matrix.....
dXdV = zeros(7,9);
dXdV(1:3,1) = 1;
dXdV(7,2) = 1;
dXdV(1:6,4:9) = eye(6);


J2sqrt = 0;
for i = 1:6
    number = 1;
    if i > 3
        number = 2;
    end
    J2sqrt = J2sqrt + number*s(i)^2;
end
J2sqrt = sqrt(0.5*J2sqrt);

dRdE = (zeros(9,6));
dummy = (-p)/(m-1)*( m -1 + ( sqrt(3)*J2sqrt/M/(-p))^m);
C1PlasticP  = m * (sqrt(3)*J2sqrt/M/(-p))^(m-1) * sqrt(3)*J2sqrt/M/(-p)^2 - dummy*(m-1)/p^2;
DeltaStrainVolPlastic = gamma*C1PlasticP;
dRdE(1,1:3) = p0*exp( -(DeltaStrainVol-DeltaStrainVolPlastic)/kappa)/kappa;
dRdE(4:9,1:6) = -2*G*(eye(6)-1/3*mIdentity*mIdentity')*ThisDev;

dRdV = Jacobian;

D = dXdV*(-dRdV\dRdE);


% now I want to compute correctly this matrix

function [Jacobian] = EvaluateJacobian(p, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev)

Jacobian = zeros(9,9);

delta = 1E-4*1i;

resRef = ComputeResidual(p, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);

res1 = ComputeResidual(p+delta, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian(:,1) = imag(res1-resRef)/imag(delta);

res1 = ComputeResidual(p, pc+delta, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian(:,2) = imag(res1-resRef)/imag(delta);

res1 = ComputeResidual(p, pc, gamma+delta, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian(:,3) = imag(res1-resRef)/imag(delta);

for i = 1:6
    s1 = s;
    s1(i) = s1(i)+delta;
    res1 = ComputeResidual(p, pc, gamma, s1, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
    Jacobian(:,3+i) = imag(res1-resRef)/imag(delta);
end


Jacobian2 = zeros(9,9);

delta = 1E-6;

resRef = ComputeResidual(p, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);

res1 = ComputeResidual(p+delta, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian2(:,1) = (res1-resRef)/(delta);

res1 = ComputeResidual(p, pc+delta, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian2(:,2) = (res1-resRef)/(delta);

res1 = ComputeResidual(p, pc, gamma+delta, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
Jacobian2(:,3) = (res1-resRef)/(delta);

for i = 1:6
    s1 = s;
    s1(i) = s1(i)+delta;
    res1 = ComputeResidual(p, pc, gamma, s1, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
    Jacobian2(:,3+i) = (res1-resRef)/(delta);
end
% Jacobian = -Jacobian2;
hola = 1;

function  [residual] = ComputeResidual(p, pc, gamma, s, p0, pc0, s0, DeltaStrainVol, DeltaStrainDev)

[kappa, lambda, M, nu, n, r, m] = GetConstitutiveParametersCASM();


residual = zeros(9,1);



ThisDev = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.0],[6,6]);



K = (-p)/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);

J2sqrt = 0;
for i = 1:6
    number = 1;
    if i > 3
        number = 2;
    end
    J2sqrt = J2sqrt + number*s(i)^2;
end
J2sqrt = sqrt(0.5*J2sqrt);

% V1 = [1;1;1;0;0;0]/3;
V2 = s;
if ( J2sqrt > 1E-5)
    V2 = V2/(2*J2sqrt);
end
V2(4:6) = 2*V2(4:6);

% C1Yield = 1/p/log(r) + n*( sqrt(3)*J2sqrt/M)^n*1/(-p)^(n+1);
% C2Yield = n* 3^(n/2)/M^n/(-p)^n * J2sqrt^(n-1);

dummy = (-p)/(m-1)*( m -1 + ( sqrt(3)*J2sqrt/M/(-p))^m);


C1PlasticP  = m * (sqrt(3)*J2sqrt/M/(-p))^(m-1) * sqrt(3)*J2sqrt/M/(-p)^2 - dummy*(m-1)/p^2;
C2PlasticP  = m * (sqrt(3)*J2sqrt/M/(-p))^(m-1) * sqrt(3)/M/(-p);

DeltaStrainVolPlastic = gamma*C1PlasticP;
DeltaStrainDevPlastic = gamma*C2PlasticP*V2;

residual(1) = p-p0*exp( -(DeltaStrainVol-DeltaStrainVolPlastic)/kappa);
residual(2) = pc - pc0 * exp( -DeltaStrainVolPlastic/(lambda-kappa));
residual(3) = EvaluateYieldSurface(p, pc, J2sqrt, M, n, r);

residual(4:9) = s - s0 - 2 * G * ThisDev * (DeltaStrainDev - DeltaStrainDevPlastic);


function yield = EvaluateYieldSurface(p, pc, J2sqrt, M, n, r)

yield = (sqrt(3)*J2sqrt/M/(-p))^n + 1/log(r)*log(p/pc);



function [D2] = CheckDerivative(X, DeltaStrain)

delta = 1E-6;

for i = 1:6
    DS = DeltaStrain;
    DS(i) = DS(i)+delta;
    [X2] = ReturnMapping(X, DS);
    
    DS = DeltaStrain;
    DS(i) = DS(i)-delta;
    [X3] = ReturnMapping(X, DS);
    
    D2(:,i) = (X2-X3)/(2*delta);
end


function [D, De, H] = ComputeContinuousMatrix(sigma, pc, plastic)

[kappa, lambda, M, nu, n, r, m] = GetConstitutiveParametersCASM();


p = mean(sigma(1:3));

J2sqrt = 0;
ss = sigma;
ss(1:3) = ss(1:3)-p;
for i = 1:6
    number = 1;
    if i > 3
        number = 2;
    end
    J2sqrt = J2sqrt + number*ss(i)^2;
end
J2sqrt = sqrt(0.5*J2sqrt);



K = (-p)/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);



De = (zeros(6,6));
De(1:3,1:3) = K;
De = De + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3); zeros(3,3), 0.5*eye(3,3)];

if (plastic == false)
    D = De; 
    H = 0;
    return;
end



dummy1 = (-p)/(m-1)*( m -1 + ( sqrt(3)*J2sqrt/M/(-p))^m);

C1F =  1/(p*log(r)) + (3^(1/2)*J2sqrt*n*(-(3^(1/2)*J2sqrt)/(M*p))^(n - 1))/(M*p^2);
C1PlasticP  = m * (sqrt(3)*J2sqrt/M/(-p))^(m-1) * sqrt(3)*J2sqrt/M/(-p)^2 - dummy1*(m-1)/p^2;

C2F = -(3^(1/2)*n*(-(3^(1/2)*J2sqrt)/(M*p))^(n - 1))/(M*p);
C2PlasticP  = m * (sqrt(3)*J2sqrt/M/(-p))^(m-1) * sqrt(3)/M/(-p);

V2 = ss;
V2(4:6) = 2*V2(4:6);
if ( J2sqrt > 1E-9)
    V2 = V2/(2*J2sqrt);
end

DeltaStrainVolPlastic = C1PlasticP;
DeltaStrainDevPlastic = C2PlasticP*V2;

gradG = DeltaStrainDevPlastic;
gradG(1:3) = gradG(1:3) + DeltaStrainVolPlastic/3;


gradF = C2F*V2;
gradF(1:3) = gradF(1:3) + C1F/3;


H = -1/(pc*log(r)) * pc/(lambda-kappa) * DeltaStrainVolPlastic;

D = De - De*gradG*gradF'*De/(H + gradF'*De*gradG);