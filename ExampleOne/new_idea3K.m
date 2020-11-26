function [] = new_idea3K()

addpath('../')
% 1. Define the problem

T = 1E-5;


CP.E = 1000;
CP.nu = 0.3;
CP.k = 1E-2;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.02;

model = createpde(1);

dx = 0.05; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';



ElementType = 'T3T3';
Nodes = Nodes;
Elements = Elements;

NSteps = 10.^[0:5];
NSteps = floor(NSteps); NSteps = sort(NSteps);





figure(99)
for i = 1:4; subplot(2,2,i); hold off; end


RKMethod = 2;

kk = 10.^[-5:1:-1];
dtdt = 10.^[-2:1:1];

[dtdt, kk] = meshgrid(dtdt, kk)
Values = nan*kk;
figure(212); hold off;

Values = nan*kk;
Values2 = nan*kk;
for auxi = 1:size(kk,1)
    
    for auxj = 1:size(kk,2)
        E = 100;
        
        k = kk(auxi, auxj);
        
        CP.E = E;
        nu = CP.nu;
        CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu)
        CP.k = k;
        
        M(auxi,auxj) = CP.M;
    end
end


for auxi = 1:size(kk,1)
    
    for auxj = 1:size(kk,2)
        E = 100;
        
        k = kk(auxi, auxj);
        
        CP.E = E;
        nu = CP.nu;
        CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
        CP.k = k;
        
        
        
        
        
        dt = 0;
        [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
        
        [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, 0);
        
        
        [C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
        
        nNodes = size(Nodes, 1);
        A = C\(K);
        dt = dtdt(auxi, auxj);
        
        ii = eye(3*nNodes, 3*nNodes);
        B = ii+dt*A;
        
        [values] = eig(full(B), 'nobalance');
        values = real(values);
        figure(212); plot(sort(values)); hold on
        Values(auxi,auxj) = max(abs(values));
        Values2(auxi,auxj) = min(abs(values));
        
        figure(22)
        contourf(log10(dtdt), log10(M.*kk), log10(Values), [0:0.5:5]);
        colorbar
        drawnow
        
        %planeFitting?
        nn = size(kk,1)*size(kk,2);
        A = [ones(nn,1) reshape( log10(M.*kk), [nn,1]), reshape( log10(dtdt), [nn,1])];
        b = reshape(log10(Values), [nn,1]);
        x = A\b
        hola = 1;
    end
end