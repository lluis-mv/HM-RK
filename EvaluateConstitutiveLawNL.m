function GPInfo = EvaluateConstitutiveLawNL(GPInfo, CP, dt, k, a, b, nn, initialize)

nElem = size(GPInfo,1);

strain= zeros(6,1);

finalize = false;

if (nargin == 6)
    initialize = false;
    finalize = true;
    nn = length(b);
end

kSigma = zeros( 7,nn);
[kappa, lambda, M, nu] = GetConstitutiveParameters();



for el = 1:nElem
    
    for gp = 1:size(GPInfo,2)
        
        sigma0 = [GPInfo(el,gp).StressPrev; GPInfo(el,gp).HistoryPrev];
        nSystem = GPInfo(el,gp).dofsU;
        
        for i = 1:nn
            sigmaStep = sigma0;
            
            for j = 1:i-1
                sigmaStep = sigmaStep + dt*a(i,j)*kSigma(:,j);
            end
            
            strain([1,2,4]) = (GPInfo(el,gp).B)*k(nSystem, i);
            if ( initialize && i == nn)
                D = ComputeD(CP, sigmaStep, strain);
            else
                D = GPInfo(el,gp).DPrev(:,:,i);
            end

            kSigma(:,i) = D*strain;
        end
        
        
        D = ComputeD(CP, sigmaStep, strain);
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
            GPInfo(el,gp).DPrev(:,:,1) = D;
            m = [1,1,1,0,0,0]';
            
%             if ( CP.Elastic)
%                 GPInfo(el, gp).ConstrainedModulus =  mean(abs(sigma(1:3)))/kappa;
%             else
%                 GPInfo(el, gp).ConstrainedModulus =  mean(abs(sigma(1:3)))/lambda * number;
%             end
            
            D = ComputeD(CP, sigma, strain);
        end
        
        GPInfo(el,gp).D6 = D(1:6,1:6);
        GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
    end
    
end



function D = ComputeD(CP, X, DeltaStrain)
X = X;
p = 1/3*(X(1)+X(2)+X(3));
[kappa, lambda, M, nu] = GetConstitutiveParameters();
K = -p/kappa;

G = 3*K*(1-2*nu)/2/(1+nu);

%De = K * [ones(3,3), zeros(3,3); zeros(3,6)] ...
%    + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3);
%    zeros(3,3), 0.5*eye(3,3)];
De = reshape([G.*(4.0./3.0)+K,G.*(-2.0./3.0)+K,G.*(-2.0./3.0)+K,0.0,0.0,0.0,G.*(-2.0./3.0)+K,G.*(4.0./3.0)+K,G.*(-2.0./3.0)+K,0.0,0.0,0.0,G.*(-2.0./3.0)+K,G.*(-2.0./3.0)+K,G.*(4.0./3.0)+K,0.0,0.0,0.0,0.0,0.0,0.0,G,0.0,0.0,0.0,0.0,0.0,0.0,G,0.0,0.0,0.0,0.0,0.0,0.0,G],[6,6]);

gradYield = GradYieldSurface(X);


loadingCondition = gradYield'*De*DeltaStrain;


if ( CP.Elastic == true)
    loadingCondition = -100;
end
% loadingCondition = 1;
if ( loadingCondition <= 0)
    D = [De; zeros(1,6)];
    return;
end




pc = X(7);

[u] = GetConstitutiveParameterU();
yield = YieldSurface(X);

R = yield/(-2*pc);

U = u/tan(pi/2*R);

tracegradYield = gradYield(1)+gradYield(2)+gradYield(3);
H = (-tracegradYield/(lambda-kappa) + U/R  ) * transpose(gradYield)*X(1:6);


D = De - De*gradYield* transpose(gradYield)*De/( H + transpose(gradYield)*De*gradYield);

D = [D;
    -pc/(lambda-kappa)* sum(gradYield(1:3))* transpose(gradYield)*De/(H + transpose(gradYield)*De*gradYield)];








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

