% Function to perform all the operations

function [Xnew, D, Dcontinuous] = ExplicitCamClay(X, DeltaStrain, method)


global tolerance
tolerance = 5E-13;

if (method > 0)
    
    [beta, dBetadEpsilon, Dinter] = GetStressIntersection(X, DeltaStrain, method);
    
    if (beta == 1)
        plastic = false;
        [Xnew, D] = IntegrateStress(X, DeltaStrain, method, false);
    elseif ( beta == 0)
        plastic = true;
        [Xnew, D] = IntegrateStress(X, DeltaStrain, method, true);
        [Xnew, D] = ReturnMapping(Xnew, DeltaStrain, false, D);
    else
        plastic = true;
        [Xnew, ~] = IntegrateStress(X, beta*DeltaStrain, method, false);
        [Xnew, D] = IntegrateStress(Xnew, (1-beta)*DeltaStrain, method, true, beta, dBetadEpsilon, Dinter);
        [Xnew, D] = ReturnMapping(Xnew, DeltaStrain, false, D);
    end
else
    [Xnew, D, plastic] = ReturnMapping(X, DeltaStrain, true);
    %     Dcontinuous = D;
    %     return;
end


[Dcontinuous] = GetContinuousConstitutiveMatrix(Xnew, plastic);


function [XNew, D, plastic] = ReturnMapping( X, DeltaStrain, returnMapping, Di)

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



if (returnMapping)
    p = p0 * exp( DeltaStrainVol/kappa);
    pc = pc0;
    
    K = p/kappa;
    G = 3*K*(1-2*nu)/2/(1+nu);
    
    s = s0 + 2*G*ThisDev * DeltaStrainDev;
else
    p = p0;
    s = s0;
    pc = pc0;
    DeltaStrainVol = 0;
    DeltaStrainDev = 0*DeltaStrainDev;
end

yield = YieldSurface( [p*mIdentity+s; pc] );

if ( returnMapping && yield < tolerance)
    XNew = [p*mIdentity+s; pc];
    
    dPdeVol = p0 * exp( DeltaStrainVol/kappa)/kappa;
    D = zeros(6,6);
    D(1:3,1:3) = dPdeVol;
    D = D + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3); zeros(3,3), 0.5*eye(3,3)];
    
    dGdEvol = 3*(1-2*nu)/2/(1+nu)*p0*exp(DeltaStrainVol/kappa)/kappa^2;
    
    D(:,1:3) = D(:,1:3)  + 2*ThisDev*DeltaStrainDev*dGdEvol;
    D = [D; zeros(1,6)];
    return;
    
elseif ( returnMapping == false && abs(yield) < tolerance)
    XNew = [p*mIdentity+s; pc];
    D = Di;
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
    
    
% % % %     if ( rcond(Jacobian) < 1E-15) || ( isnan(rcond(Jacobian) ))
% % % %         %         [Jacobian, p, pc, s, gamma] = RMLineSearch(kappa, lambda, nu, M, ...
% % % %         %          p0, pc0, s0, DeltaStrainVol, DeltaStrainDev);
% % % %         hola = 1;
% % % %     end
    
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
K = p/kappa;
G = 3*K*(1-2*nu)/2/(1+nu);

XNew = [p*mIdentity+s; pc];


% Now we have to obtain the stiffness matrix.....
dXdV = zeros(7,9);
dXdV(1:3,1) = 1;
dXdV(7,2) = 1;
dXdV(1:6,4:9) = eye(6);

dVdX = zeros(9,7);
dVdX(1,1:3) = 1/3;
dVdX(2,7) = 1;
Idev = eye(6);
Idev(1:3,1:3) = Idev(1:3,1:3)-1/3;
dVdX(4:9,1:6) = Idev;

dRdV0 = (zeros(9,9));
dRdV0(1,1) = -exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa);
dRdV0(2,2) = -exp( 2 * gamma * (p-pc) /(lambda-kappa));
dRdV0(4:9,4:9) = -eye(6);


dRdE = (zeros(9,6));
dRdE(1,1:3) = -p0*exp( (DeltaStrainVol-2*gamma*(p-pc) )/kappa)/kappa;
dRdE(4:9,1) =  6*gamma/M^2* 3*(1-2*nu)/2/(1+nu)/kappa;
dRdE(4:9,1:6) = -2*G*(eye(6)-1/3*mIdentity*mIdentity')*ThisDev;

dRdV = Jacobian;

if ( returnMapping == true)
    Di = zeros(7,6);
else
    dRdE = 0*dRdE;
end

dVdE = -(dRdV)\(  dRdV0*dVdX*Di + dRdE);
D = dXdV*dVdE;


function [beta, dBetadEpsilon, D] = GetStressIntersection(X, DeltaStrain, method)

global tolerance
dBetadEpsilon = [];
D = [];

yield0 = YieldSurface(X);

Xend = IntegrateStress(X, DeltaStrain, method, false);
yieldEnd = YieldSurface(Xend);

if ( yieldEnd < tolerance)
    beta = 1;
    return;
end



if ( yield0 > -tolerance)
    
    dSigma = Xend-X; dSigma = dSigma(1:6);
    gradF = GradYieldSurface(X);
    if ( dSigma'*gradF > 0 || norm(DeltaStrain) == 0 )
        beta = 0;
        return;
    end
    
end


% now with a good method
beta = 1;

normE = 100;

dBeta = 1;
iter = 0;
while ( normE > tolerance)
    
    [Xbeta, D] = IntegrateStress(X, beta*DeltaStrain, method, false);
    yieldbeta = YieldSurface(Xbeta);
    
    normE = norm(yieldbeta);
    if ( normE < tolerance && abs(dBeta) < tolerance)
        break;
    elseif ( abs(dBeta) < tolerance && normE < 100*tolerance)
        break;
    elseif ( iter == 25)
        betas = 0.05:0.1:0.95;
        for i = 1:length(betas)
            [Xbeta, D] = IntegrateStress(X, betas(i)*DeltaStrain, method, false);
            yieldbeta(i) = YieldSurface(Xbeta);
            normE(i) = norm(yieldbeta(i));
            if (yieldbeta(i) > tolerance)
                beta = betas(i);
                break;
            end
        end
        iter = iter+1;
        continue;
    elseif ( iter == 50)
        betas = 0.01:0.025:0.975;
        for i = 1:length(betas)
            [Xbeta, D] = IntegrateStress(X, betas(i)*DeltaStrain, method, false);
            yieldbeta(i) = YieldSurface(Xbeta);
            normE(i) = norm(yieldbeta(i));
            if (yieldbeta(i) > tolerance)
                beta = betas(i);
                break;
            end
        end
        iter = iter+1;
        continue;
    elseif ( iter > 100)
        disp(['The find beta is not working: ', num2str(normE)]);
        
        beta = 1;
        break
    end
    
    
    gradF = GradYieldSurface(Xbeta);
    
    gradient = gradF(1:6)'*D(1:6,1:6)*DeltaStrain;
    
    dBeta = -yieldbeta/gradient;
    beta = beta + dBeta;
    if ( beta > 1.01)
        beta = 0.25;
    elseif (beta < -0.01)
        beta = 0.75;
    end
    iter = iter+1;
end

[Xbeta, ~, A] = IntegrateStress(X, beta*DeltaStrain, method, false);
gradF = GradYieldSurface(Xbeta);



dBetadEpsilon = zeros(6,1);
dKdET = zeros(6,6,length(A.bRK) );
ind = 1:6;
iter = 0;

while ( iter < 100)
    for i = 1:length(A.bRK)
        term = zeros(6,6);
        for j = 1:i-1
            term = term + A.aRK(i,j) * dKdET(ind,ind,j);
        end
        dKdET(:,:,i) = A.dKdX(ind,ind,i)*term + A.dKdE(ind,ind,i)*( beta*eye(6) + DeltaStrain*dBetadEpsilon');
    end
    
    denom = 0;
    for i = 1:length(A.bRK)
        denom = denom + A.bRK(i)*gradF'* A.dKdE(ind,ind,i)*DeltaStrain;
    end
    
    numerator = 0*gradF;
    
    for i = 1:length(A.bRK)
        term = zeros(6,6);
        for j = 1:i-1
            term = term + A.aRK(i,j) * dKdET(ind,ind,j);
        end
        aux =  A.dKdX(ind,ind,i)*term + A.dKdE(ind,ind,i)*beta;
        numerator = numerator + A.bRK(i)* (gradF'*aux)';
    end
    
    dBetadEpsilonNew = - numerator/denom;
    normdA = norm(dBetadEpsilonNew-dBetadEpsilon);
    if ( normdA < tolerance)
        break;
    end
    if ( iter > 200)
        error('Not capable of finding the derivative')
    end
    dBetadEpsilon = dBetadEpsilonNew;
    iter = iter+1;
end



D = zeros(6,6);
for i = 1:length(A.bRK)
    term = zeros(6,6);
    for j = 1:i-1
        term = term + A.aRK(i,j) * dKdET(ind,ind,j);
    end
    dKdET(:,:,i) = A.dKdX(ind,ind,i)*term + A.dKdE(ind,ind,i)*( beta*eye(6) + DeltaStrain*dBetadEpsilon');
    D = D + A.bRK(i)*dKdET(:,:,i);
end


D = [D; zeros(1,6)];

function [Xnew, D, AllTerms] = IntegrateStress(X, DeltaStrain, method, plastic, beta, dBetadEpsilon, Dinter)

[aRK, bRK] = GetRungeKutta(method);

k = zeros( length(X), length(bRK));

dKdX = zeros(7,7,length(bRK));
dKdE = zeros(7,6,length(bRK));

for i = 1:length(bRK)
    Xstep = X;
    for j = 1:i-1
        Xstep = Xstep + aRK(i,j)*k(:,j);
    end
    
    [k(:,i)] = SourceTerm(Xstep, DeltaStrain, plastic);
    [ dKdX(:,:,i), dKdE(:,:,i)] = SourceTermDerivatives( Xstep, DeltaStrain, plastic);
    
end

Xnew = X;
for i = 1:length(bRK)
    Xnew = Xnew + bRK(i)*k(:,i);
end

if ( nargin == 4)
    dKdET = zeros(7,6,length(bRK));
    D = zeros(7,6);
    for i = 1:length(bRK)
        term = zeros(7,6);
        for j = 1:i-1
            term = term + aRK(i,j)*dKdET(:,:,j);
        end
        dKdET(:,:,i) = dKdX(:,:,i)*term + dKdE(:,:,i);
        D = D + bRK(i)*dKdET(:,:,i);
    end
else
    
    % The same but with more terms....
    dKdET = zeros(7,6,length(bRK));
    thisI = eye(6,6);
    D = zeros(7,6);
    for i = 1:length(bRK)
        term = zeros(7,6);
        for j = 1:i-1
            term = term + aRK(i,j)*dKdET(:,:,j);
        end
        dKdET(:,:,i) = dKdX(:,:,i)*(Dinter + term) + dKdE(:,:,i)* ( (1-beta)*thisI + DeltaStrain*(-dBetadEpsilon')/(1-beta) ) ;
        D = D + bRK(i)*dKdET(:,:,i);
    end
    D = D + Dinter;
    
    
end
if (nargout == 3)
    AllTerms.bRK = bRK;
    AllTerms.aRK = aRK;
    AllTerms.dKdX = dKdX;
    AllTerms.dKdE = dKdE;
end


function [S] = SourceTerm(X, DeltaStrain, plastic)

[kappa, lambda] = GetConstitutiveParameters();

S = zeros(7,1);

sigma = X(1:6);
% p = mean(sigma(1:3));
p = 1/3*(sigma(1)+sigma(2)+sigma(3));
pc = X(7);
De = GetDe(p);

if ( plastic == false)
    S(1:6) = De*DeltaStrain;
    return;
end


H = 2*p * pc/(lambda-kappa) * 2 * (p-pc);
gradYield = GradYieldSurface(X);

Dep = De - De*gradYield* transpose(gradYield)*De/( H + transpose(gradYield)*De*gradYield);
S(1:6) = Dep*DeltaStrain;

IncrGamma = transpose(gradYield)*De*DeltaStrain/(H + transpose(gradYield)*De*gradYield);
S(7) = pc/(lambda-kappa)* sum(gradYield(1:3))*IncrGamma;



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




function [partialKpartialX, partialKpartialE] = SourceTermDerivatives(X, DeltaStrain, plastic)

global NumericalDerivative;
if (isempty(NumericalDerivative) )
    NumericalDerivative = true;
end

if (NumericalDerivative == false)
    [partialKpartialX, partialKpartialE] = GetPartialDerivatives(X, DeltaStrain,plastic);
    return
end


partialKpartialX = zeros(7,7);
partialKpartialE = zeros(7,6);



delta = 1E-6*1i;

for i = 1:7
    X1 = X;
    X1(i) = X1(i)+delta;
    [S1] = SourceTerm(X1, DeltaStrain, plastic);
    partialKpartialX(:,i) = imag(S1)/imag(delta);
end



for i = 1:6
    DeltaStrain1 = DeltaStrain;
    DeltaStrain1(i) = DeltaStrain1(i)+delta;
    [S1] = SourceTerm(X, DeltaStrain1, plastic);
    partialKpartialE(:,i) = imag(S1)/imag(delta);
end



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

dYield = 2/3*(p-pc)*[m] + 2*3/M^2 * jj22 * V2;
