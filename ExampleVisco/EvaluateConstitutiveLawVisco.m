function GPInfo = EvaluateConstitutiveLawVisco(GPInfo, CP, dt, k, a, b, nn, initialize)

nElem = size(GPInfo,1);

strain= zeros(6,1);

finalize = false;

if (nargin == 6)
    initialize = false;
    finalize = true;
    nn = length(b);
end

kSigma = zeros( 6,nn);
[kappa, lambda, M, nu] = GetConstitutiveParameters();

nVisco = 1;
mVisco = 0.02;

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
                eVP = zeros(6,1); 
                
                if ( YieldSurface(sigmaStep) > 0)
                    eVP = GradYieldSurface(sigmaStep);
                    eVP = eVP*(YieldSurface(sigmaStep)/mVisco)^nVisco;
                end
                
            else
                D = GPInfo(el,gp).DPrev(:,:,i);
                eVP = GPInfo(el,gp).eVPPrev(:,i);
            end

            kSigma(:,i) = D*(strain-eVP);
        end
        
        
        D = ComputeD(CP, sigmaStep, strain);
        GPInfo(el,gp).DPrev(:,:,nn) = D;
        if ( YieldSurface(sigmaStep) > 0)
            eVP = GradYieldSurface(sigmaStep);
            eVP = eVP*(YieldSurface(sigmaStep)/mVisco)^nVisco;
        end
        GPInfo(el,gp).eVPPrev(:,nn) = eVP;
        
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
            GPInfo(el,gp).eVPPrev(:,1) = eVP;
            m = [1,1,1,0,0,0]';
            
            
            D = ComputeD(CP, sigma, strain);
        end
        
        GPInfo(el,gp).D6 = D(1:6,1:6);
        GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
        GPInfo(el,gp).eVP = eVP;
    end
    
end



function D = ComputeD(CP, X, DeltaStrain)
X = X;
p = 1/3*(X(1)+X(2)+X(3));
[kappa, lambda, M, nu] = GetConstitutiveParameters();

E = 4000; 
nu = 0.3;
K = E/3/(1-2*nu); 
G = E/2/(1+nu);


D = reshape([G.*(4.0./3.0)+K,G.*(-2.0./3.0)+K,G.*(-2.0./3.0)+K,0.0,0.0,0.0,G.*(-2.0./3.0)+K,G.*(4.0./3.0)+K,G.*(-2.0./3.0)+K,0.0,0.0,0.0,G.*(-2.0./3.0)+K,G.*(-2.0./3.0)+K,G.*(4.0./3.0)+K,0.0,0.0,0.0,0.0,0.0,0.0,G,0.0,0.0,0.0,0.0,0.0,0.0,G,0.0,0.0,0.0,0.0,0.0,0.0,G],[6,6]);



function yield = YieldSurface(X)

Su = 10;

yield = 0;
sigma = X(1:6);


p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    yield = yield + ( sigma(i)-p)^2;
end
for i = 4:6
    yield = yield + 2*(sigma(i))^2;
end
yield = sqrt(3)*sqrt(0.5*yield)-2*Su;
yield = yield/(2*Su);

function dYield = GradYieldSurface(X)
Su = 10;

jj22 = 0;
sigma = X(1:6);

% p = mean(sigma(1:3));
p = 1/3*(sigma(1)+sigma(2)+sigma(3));
for i = 1:3
    jj22 = jj22 + ( sigma(i)-p)^2;
end
for i = 4:6
    jj22 = jj22 + 2*(sigma(i))^2;
end
jj22 = sqrt(0.5*jj22);

V2 = sigma; V2(1:3) = V2(1:3)-p*ones(3,1); V2(4:6) = 2*V2(4:6);
if ( norm(jj22) > 1E-12)
    V2 = V2/(2*jj22);
else
    V2 = zeros(6,1);
end

dYield = sqrt(3)*V2;
    
dYield = dYield/(2*Su);
