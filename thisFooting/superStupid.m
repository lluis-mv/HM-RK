function [] = superStupid


X0 = [10; 20];

RKRef = 8;
RKMethods = [8, 1:7];

NSTEPS = 2.^(0:1:12);

index = 1;


% figure1
X = X0;
nSteps = 100;
x = zeros(1,nSteps);
x(1) = X(1);
for i = 1:nSteps
    
    X = ComputeThis(X, 1, 8, 1/nSteps);
    x(i+1) = X(1);
end
figure(4)
semilogx(x, -(0:nSteps)/nSteps, '-.')

for nSteps = NSTEPS
    for RK = RKMethods    
        X = ComputeThis(X0, nSteps, RK);
        N(index, RK) = X;
        if ( RK == RKRef)
            Nadim = X;
        end
    end
    
    figure(1)
    clf
    for i = 1:size(N,2)
        loglog( 1./NSTEPS(1:index), N(1:index,i), '*-.')
        hold on
    end
    grid
    
    figure(2)
    clf
    for i = 1:size(N,2)
        loglog( 1./NSTEPS(1:index), abs(N(1:index,i)-Nadim), '*-.')
        hold on
    end
    grid
    index = index+1;
end


function [X] = ComputeThis(X, nSteps, RK, DeltaStrain)

if (nargin == 3)
    DeltaStrain = 1;
end
[a,b] = GetRungeKutta(RK);
k = zeros(  2, length(b));

for step = 1:nSteps
    
    for i = 1:length(b)
        XStep = X;
        for j = 1:i-1
            XStep = XStep + a(i,j)*k(:,j);
        end
        
        k(:,i) = SourceTerm(XStep)*DeltaStrain/nSteps;
    end
    
    for i = 1:length(b)
        X = X + b(i)*k(:,i);
    end
end
if (nargin == 3)
    X = X(1);
end

function [y] = SourceTerm(X)
kappa = 0.05;
lambda = 0.5;
p = X(1);
pc = X(2);
R = p/(2*pc);
u = 10;
U = u/tan(pi/2*R);

De = p/kappa;

H = (1/(lambda-kappa) + U/R)*p;
Dep = De - De*1*1*De/(H+De);
Dh = pc/(lambda-kappa)*De/(H+De);
y = [Dep;Dh];



