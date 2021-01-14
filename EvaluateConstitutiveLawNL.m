function GPInfo = EvaluateConstitutiveLawNL(GPInfo, X, dt, k, a, b, c, nn, initialize)

nElem = size(GPInfo,1);

strain= zeros(6,1);

finalize = true;
if (nargin >= 8)
    finalize = false;
else
    nn = length(b);
    initialize = false;
end
if (nargin < 9)
    initialize = false;
end


kSigma = zeros( 7,nn);


% if (nargin == 8)
for el = 1:nElem
    for gp = 1:size(GPInfo,2)
        
        if ( el == 1 && gp == 1 && initialize == false)
            hola = 1;
        end
        sigma0 = [GPInfo(el,gp).StressPrev; GPInfo(el,gp).HistoryPrev];
        nSystem = GPInfo(el,gp).dofsU;
        for i = 1:nn
            sigmaStep = sigma0;
            
            for j = 1:i-1
                sigmaStep = sigmaStep + dt*a(i,j)*kSigma(:,j);
            end
            
            strain([1,2,4]) = (GPInfo(el,gp).B)*k(nSystem, i);
            D = ComputeD(sigmaStep, strain);
            if ( initialize == true)
                D = ComputeD(sigmaStep, strain);
            elseif ( finalize == false && i < nn)
                D = GPInfo(el,gp).DPrev(:,:,i);
            elseif ( finalize)
                D = GPInfo(el,gp).DPrev(:,:,i);
            end
            kSigma(:,i) = D*strain;
        end
        D = ComputeD(sigmaStep, strain);
        GPInfo(el,gp).DPrev(:,:,nn) = D;
        
        if ( finalize)
            sigma = sigma0;
            for i = 1:length(b)
                sigma = sigma + dt*b(i)*kSigma(:,i);
            end
            GPInfo(el,gp).StressNew = sigma(1:6);
            GPInfo(el,gp).StressPrev = sigma(1:6);
            GPInfo(el,gp).HistoryNew = sigma(7:end);
            GPInfo(el,gp).HistoryPrev = sigma(7:end);
            D = ComputeD(sigma, strain);
            GPInfo(el,gp).DFin = D;
            GPInfo(el,gp).DPrev(:,:,1) = D;
            m = [1,1,1,0,0,0]';
            GPInfo(i, gp).ConstrainedModulus =  m'*D(1:6,1:6)*m/100;
        end
        if ( initialize)
            GPInfo(el,gp).DFin = D;
            GPInfo(el,gp).DPrev(:,:,nn) = D;
        end
        
        GPInfo(el,gp).D6 = D(1:6,1:6);
        GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
    end
    
end

% else
%     for el = 1:nElem
%         for gp = 1:size(GPInfo,2)
%
%
%             sigma0 = GPInfo(el,gp).StressPrev;
%             nSystem = GPInfo(el,gp).dofsU;
%
%             for i = 1:length(b)
%                 sigmaStep = sigma0;
%
%                 for j = 1:i-1
%                     sigmaStep = sigmaStep + dt*a(i,j)*kSigma(:,j);
%                 end
%                 D = ComputeD(sigmaStep);
%                 strain([1,2,4]) = (GPInfo(el,gp).B)*k(nSystem, i);
%                 kSigma(:,i) = D*strain;
%             end
%
%
%             GPInfo(el,gp).D6 = D;
%             GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
%         end
%
%     end
% end


function De = ComputeD(X, DeltaStrain)
X = X;
p = 1/3*(X(1)+X(2)+X(3));
[kappa, lambda, M, nu] = GetConstitutiveParameters();
K = -p/kappa;

G = 3*K*(1-2*nu)/2/(1+nu);

De = K * [ones(3,3), zeros(3,3); zeros(3,6)] ...
    + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3);
    zeros(3,3), 0.5*eye(3,3)];

gradYield = GradYieldSurface(X);


loadingCondition = gradYield'*De*DeltaStrain;




if ( loadingCondition < 0)
    De = [De; zeros(1,6)];
    return;
end




pc = X(7);

[u] = GetConstitutiveParameterU();
yield = YieldSurface(X);

R = yield/(-2*pc);

U = u/tan(pi/2*R);

tracegradYield = gradYield(1)+gradYield(2)+gradYield(3);
H = (-tracegradYield/(lambda-kappa) + U/R  ) * transpose(gradYield)*X(1:6);


De = De - De*gradYield* transpose(gradYield)*De/( H + transpose(gradYield)*De*gradYield);

De = [De;
    -pc/(lambda-kappa)* sum(gradYield(1:3))* transpose(gradYield)*De/(H + transpose(gradYield)*De*gradYield)];
hola = 1;








function yield = YieldSurface(X)
[~, ~, M, ~] = GetConstitutiveParameters();
jj22 = 0;
sigma = X(1:6);

p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    jj22 = jj22 + ( sigma(i)-p)^2;
end
for i = 4:6
    jj22 = jj22 + 2*(sigma(i))^2;
end
p = -p;
jj22 = sqrt(0.5*jj22);
yield = p*((3*jj22^2)/(M^2*p^2) + 1);



function dYield = GradYieldSurface(X)
[~, ~, M, ~] = GetConstitutiveParameters();
jj22 = 0;
sigma = X(1:6);
pc = X(7);
m = [1;1;1;0;0;0];


p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    jj22 = jj22 + ( sigma(i)-p)^2;
end
for i = 4:6
    jj22 = jj22 + 2*(sigma(i))^2;
end
jj22 = sqrt(0.5*jj22);


V2 = (sigma-p*m);
for i = 4:6
    V2(i) = 2*V2(i);
end

if ( norm(jj22) > 1E-12)
    V2 = V2/(2*jj22);
else
    V2 = zeros(6,1);
end

C1 = 1 - (3*jj22^2)/(M^2*p^2);
C2 = -(6*jj22)/(M^2*p);
dYield = -C1*m/3 + C2 * V2;

normJ = 0;
for i = 1:3
    normJ = normJ + dYield(i)^2;
end

for i = 4:6
    normJ = normJ + 2*(dYield(i)/2)^2;
end
dYield = dYield / normJ;

