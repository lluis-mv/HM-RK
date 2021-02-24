
function [] = main()


addpath('../')
RKMethods = [1:8];

adim = 1;
figure(99)
epsilon = 0.3;
X = 1;
for RK = RKMethods
    i = 1;
    for nSteps = 2.^[1:8]
        XNew = IntStress( 10*nSteps+ RK, X, epsilon);
        RR(i,RK) = XNew;
        if ( RK == 8)
            adim=RR(i,RK);
        end
        ddtt(i,RK) = 1/nSteps;
        loglog( ddtt, abs(RR-adim), '*-.')
        i = i+1;
    end
end


adim = 1;
figure(1)
for RK =RKMethods
    i = 1;
    for nSteps = 2.^[1:8]
        dt = 1/nSteps;
        X = ComputeThisProblem(RK, nSteps, dt);
        RR(i,RK) = X(2);
        if ( RK == 8)
            adim=RR(i,RK);
        end
        ddtt(i,RK) = dt;
        loglog( ddtt, abs(RR-adim), '*-.')
        i = i+1;
    end
end

figure(2)
% adim = 1;
for RK = RKMethods
    i = 1;
    for nSteps = 2.^[1:8]
        dt = 1/nSteps;
        X = ComputeThisProblem2(RK, nSteps, dt);
        RR(i,RK) = X(2);
        ddtt(i,RK) = dt;
        loglog( ddtt, abs(RR-adim), '*-.')
        i = i+1;
    end
end


function [X] = ComputeThisProblem(RKMethod, nSteps, dt)

X = [0;0;1];


[a,b,c] = GetRungeKutta(RKMethod);

k = zeros(  3, length(b));

B = [-1; 1];

for loadStep = 1:nSteps
    
    for i = 1:length(b)
        XStep = X;
        t = (loadStep-1)*dt + c(i)*dt;
        f = ComputeForce(t);
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
        
        D = SiffnessMatrix(XStep(3));
        C = [ zeros(2,2), B;
            D*B',    -1];
        
        C(1,:) = 0;
        C(1,1) = 1;
        
        k(:,i) = C\(f);
    end
    XNew  = X;
    for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end
    X = XNew;
end




function [X] = ComputeThisProblem2(RKMethod, nSteps, dt)

X = [0;0];


[a,b,c] = GetRungeKutta(RKMethod);

k = zeros(  2, length(b));

B = [-1; 1];

sigma = 1;

RKStress = RKMethod;
for loadStep = 1:nSteps
    
    for i = 1:length(b)
        XStep = X;
        t = (loadStep-1)*dt + c(i)*dt;
        f = ComputeForce(t);
        f = f(1:2);
        for j = 1:i-1
            XStep = XStep + dt*a(i,j)*k(:,j);
        end
        
        sigmaNew = IntStress( RKStress, sigma, B'*(XStep-X));
        D = SiffnessMatrix(sigmaNew);
        C = [B*D*B'];
        
        C(1,:) = 0;
        C(1,1) = 1;
        
        k(:,i) = C\(f);
    end
    
    XNew  = X;
    for i = 1:length(b)
        XNew = XNew + dt*b(i)*k(:,i);
    end
    sigmaNew = IntStress( RKStress, sigma, B'*(XNew-X));
    X = XNew;
    
    sigma = sigmaNew;
end

X = [X; sigma];


function XNew = IntStress( RKMethod, X, epsilon)
[a,b] = GetRungeKutta(RKMethod);
k = zeros(  1, length(b));
for i = 1:length(b)
    XStep = X;
    for j = 1:i-1
        XStep = XStep + a(i,j)*k(:,j);
    end
    
    D = SiffnessMatrix(XStep);
    
    
    k(:,i) = D*epsilon;
end
XNew = X;
for i = 1:length(b)
    XNew = XNew + b(i)*k(:,i);
end

% XNew = X*exp(epsilon/0.05);
% 
% if ( XNew > 10)
%     XNew = X * exp(epsilon/0.1);
% end

function [D] = SiffnessMatrix(sigma)
D = sigma/0.05;
if ( sigma > 10)
    D = sigma/0.1;
end


function f = ComputeForce(t)
f = [0;20*t^8;0];





