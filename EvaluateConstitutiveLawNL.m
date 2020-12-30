function GPInfo = EvaluateConstitutiveLawNL(GPInfo, X, dt, k, a, b, c, nn)

nElem = size(GPInfo,1);

strain= zeros(6,1);
if (nargin == 8)
    
    for el = 1:nElem
        for gp = 1:size(GPInfo,2)
            sigma0 = GPInfo(el,gp).StressPrev;
            nSystem = GPInfo(el,gp).dofsU;
            for i = 1:nn
                sigmaStep = sigma0;
                
                for j = 1:i-1
                    sigmaStep = sigmaStep + dt*a(i,j)*kSigma(:,j);
                end
                D = ComputeD(sigmaStep);
                strain([1,2,4]) = (GPInfo(el,gp).B)*k(nSystem, i);
                kSigma(:,i) = D*strain;
            end
            GPInfo(el,gp).D6 = D;
            GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
        end
        
    end
    
else
    for el = 1:nElem
        for gp = 1:size(GPInfo,2)
            
            
            sigma0 = GPInfo(el,gp).StressPrev;
            nSystem = GPInfo(el,gp).dofsU;
            
            for i = 1:length(b)
                sigmaStep = sigma0;
                
                for j = 1:i-1
                    sigmaStep = sigmaStep + dt*a(i,j)*kSigma(:,j);
                end
                D = ComputeD(sigmaStep);
                strain([1,2,4]) = (GPInfo(el,gp).B)*k(nSystem, i);
                kSigma(:,i) = D*strain;
            end
            
            sigma = sigma0;
            for i = 1:length(b)
                sigma = sigma + dt*b(i)*kSigma(:,i);
            end
            GPInfo(el,gp).StressNew = sigma;
            GPInfo(el,gp).StressPrev = sigma;
            D = ComputeD(sigma);
            GPInfo(el,gp).D6 = D;
            GPInfo(el,gp).D = D([1,2,4], [1,2,4]);
        end
        
    end
end


function De = ComputeD(sigmaStep)
p = 1/3*(sigmaStep(1)+sigmaStep(2)+sigmaStep(3));
[kappa, ~, ~, nu] = GetConstitutiveParameters();
K = -p/kappa;

G = 3*K*(1-2*nu)/2/(1+nu);

De = K * [ones(3,3), zeros(3,3); zeros(3,6)] ...
    + 2*G*[eye(3,3)-1/3*ones(3,3), zeros(3,3);
    zeros(3,3), 0.5*eye(3,3)];