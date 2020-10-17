% Functions to perform the implicit integration of the Modified Cam Clay

function [Xnew, Dcontinuous] = ImplicitCamClay(X, DeltaStrain)


global tolerance
tolerance = 5E-13;


[Xnew, plastic] = ReturnMapping(X, DeltaStrain);


[Dcontinuous] = GetContinuousConstitutiveMatrix(Xnew, plastic);


function [XNew, plastic] = ReturnMapping( X, DeltaStrain)

global tolerance

plastic = false;

[kappa, lambda, M, nu] = GetConstitutiveParameters();

mIdentity = [ones(3,1); zeros(3,1)];


ThisDev = eye(6); ThisDev(4:6,4:6) = ThisDev(4:6,4:6)/2;

DeltaStrainVol = sum(DeltaStrain(1:3));
DeltaStrainDev = DeltaStrain - DeltaStrainVol/3*mIdentity;

p0 = sum(X(1:3))/3;
s0 = X(1:6)-p0*mIdentity;
pc0 = X(7);




p = p0 * exp( DeltaStrainVol/kappa);
pc = pc0;

K = p/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);

s = s0 + 2*G*ThisDev * DeltaStrainDev;


yield = YieldSurface( [p*mIdentity+s; pc] );

if ( yield < tolerance)
    XNew = [p*mIdentity+s; pc];
    return;
end


plastic = true;


gamma = 0;
iter = 0;
residual = zeros(9,1);
Jacobian = zeros(9,9);
Idev =  eye(6) - mIdentity*mIdentity'/3 ;
Idev(4:6,4:6) = 2*Idev(4:6,4:6);
while (true)
    
    K = p/kappa;
    G = 3*K*(1-2*nu)/2/(1+nu);
    
    residual(1) = p-p0*exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa);
    residual(2) = pc - pc0 * exp( 2 * gamma * (p-pc) /(lambda-kappa));
    residual(3) = 3/2/M^2 * s' * Idev*s + (p-pc)^2 - pc^2;
    residual(4:9) = ( 1 + 6*G*gamma/M^2)*s - s0 - 2 * G * ThisDev * DeltaStrainDev;
    
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
    
    
    Jacobian(1,1) = 1 + 2 * gamma/kappa* p0*exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa);
    Jacobian(1,2) = - 2*gamma/kappa * p0*exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa);
    Jacobian(1,3) = 2*(p-pc)/kappa  * p0*exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa);
    
    Jacobian(2,1) = - 2*gamma/(lambda-kappa) * pc0 * exp( 2 * gamma * (p-pc) /(lambda-kappa) );
    Jacobian(2,2) = 1 + 2 * gamma/(lambda-kappa) * pc0 * exp( 2 * gamma * (p-pc) /(lambda-kappa) ) ;
    Jacobian(2,3) = -2*(p-pc)/(lambda-kappa)* pc0 * exp( 2 * gamma * (p-pc) /(lambda-kappa) ) ;
    
    Jacobian(3,1) = 2*(p-pc);
    Jacobian(3,2) = -2*p;
    thisV = s; thisV(4:6) = 2*thisV(4:6);
    Jacobian(3,4:9) = 3/M^2*thisV;
    
    
    Jacobian(4:6,1) = (6*gamma/M^2* s(1:3) - 2 * DeltaStrainDev(1:3)) * 3/kappa*(1-2*nu)/2/(1+nu) ;
    Jacobian(7:9,1) = (6*gamma/M^2* s(4:6) -  DeltaStrainDev(4:6)) * 3/kappa*(1-2*nu)/2/(1+nu) ;
    
    Jacobian(4:9,3) = 6*G/M^2*s;
    Jacobian(4:9,4:9) = (1+6*G*gamma/M^2)*eye(6);
    
    dX = -Jacobian\residual;
    
    
    p = p + dX(1);
    pc = pc+dX(2);
    gamma = gamma+dX(3);
    s = s + dX(4:9);
    
    if ( any(isnan(dX)))
        break;
    end
    iter = iter+1;
    
end

XNew = [p*mIdentity+s; pc];






function [D] = GetContinuousConstitutiveMatrix(X, plastic)

[kappa, lambda] = GetConstitutiveParameters();



sigma = X(1:6);
% p = mean(sigma(1:3));
p = 1/3*(sigma(1)+sigma(2)+sigma(3));
pc = X(7);
De = GetDe(p);

if ( plastic == false)
    D = De;
    return;
end


H = 2*p * pc/(lambda-kappa) * 2 * (p-pc);
gradYield = GradYieldSurface(X);

D = De - De*gradYield* transpose(gradYield)*De/( H + transpose(gradYield)*De*gradYield);




function De = GetDe(p)

[kappa, ~, ~, nu] = GetConstitutiveParameters();
% K = p/kappa;
% G = 3*K*(1-2*nu)/2/(1+nu);
%
% De = K * [ones(3,3), zeros(3,3); zeros(3,6)] ...
%     + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3);
%     zeros(3,3), 0.5*eye(3,3)];

De = reshape([p./kappa-(p.*(nu.*2.0-1.0).*2.0)./(kappa.*(nu+1.0)),p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),0.0,0.0,0.0,p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),p./kappa-(p.*(nu.*2.0-1.0).*2.0)./(kappa.*(nu+1.0)),p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),0.0,0.0,0.0,p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),p./kappa+(p.*(nu.*2.0-1.0))./(kappa.*(nu+1.0)),p./kappa-(p.*(nu.*2.0-1.0).*2.0)./(kappa.*(nu+1.0)),0.0,0.0,0.0,0.0,0.0,0.0,(p.*(nu.*2.0-1.0).*(-3.0./2.0))./(kappa.*(nu+1.0)),0.0,0.0,0.0,0.0,0.0,0.0,(p.*(nu.*2.0-1.0).*(-3.0./2.0))./(kappa.*(nu+1.0)),0.0,0.0,0.0,0.0,0.0,0.0,(p.*(nu.*2.0-1.0).*(-3.0./2.0))./(kappa.*(nu+1.0))],[6,6]);



function yield = YieldSurface(X)
[kappa, lambda, M, nu] = GetConstitutiveParameters();
yield = 0;
sigma = X(1:6);
pc = X(7);
% p = mean(sigma(1:3));
p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    yield = yield + ( sigma(i)-p)^2;
end
for i = 4:6
    yield = yield + 2*(sigma(i))^2;
end
yield = sqrt(0.5*yield);
yield = (sqrt(3)*yield/M)^2 + (p-pc)^2-pc^2;


function dYield = GradYieldSurface(X)
[~, ~, M, ~] = GetConstitutiveParameters();
jj22 = 0;
sigma = X(1:6);
pc = X(7);
m = [1;1;1;0;0;0];

% p = mean(sigma(1:3));
p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    jj22 = jj22 + ( sigma(i)-p)^2;
end
for i = 4:6
    jj22 = jj22 + 2*(sigma(i))^2;
end
jj22 = sqrt(0.5*jj22);

%V2 = sigma; V2(1:3) = V2(1:3)-p*ones(3,1); V2(4:6) = 2*V2(4:6);
V2 = (sigma-p*m);
for i = 4:6; V2(i) = 2*V2(i); end

if ( norm(jj22) > 1E-12)
    V2 = V2/(2*jj22);
else
    V2 = zeros(6,1);
end

dYield = 2/3*(p-pc)*m + 2*3/M^2 * jj22 * V2;
