function [] = new_idea3K()

addpath('../')
% 1. Define the problem


kk = 10.^[-3:0.5:-1]*100;
dtdt = 10.^[-4:1:1];
EE = 10.^[3:0.5:6];
hh = [0.01, 0.012, 0.015, 0.018];
hh  = [0.01, 0.015, 0.02, 0.025];



hh = sort(hh, 'descend');
dx = 0.07; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';

ElementType = 'T3T3';
RKMethod = 1;
Kc = [];
Dtc = [];
Mc = [];
Hec = [];
Valc = [];
for l = 1:length(hh)
    for i = 1:length(kk)
        for j = 1:length(dtdt)
            for ki = 1:length(EE)
                tic
                k = kk(i);
                dt = dtdt(j);
                E = EE(ki);
                he = hh(l);
                
                
                
                model = createpde(1);
                g = decsg(R1);
                geometryFromEdges(model, g);
                
                if ( all(ElementType == 'T3T3') )
                    mesh = generateMesh(model, 'Hmax', he, 'GeometricOrder','linear');
                else
                    mesh = generateMesh(model, 'Hmax', he);
                end
                
                Nodes = mesh.Nodes';
                Elements = mesh.Elements';
                
                
                
                
                
                
                CP.E = E;
                CP.nu = 0.3;
                nu = CP.nu;
                CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
                CP.k = k;
                M = CP.M;
                
                [GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
                if ( all(ElementType == 'T3T3') )
                    he = mean( sqrt([GPInfo.Weight]));
                else
                    he = 3*mean(([GPInfo.Weight]));
                end
                    
                
                [C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, 0);
                
                
                [C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
                
                nNodes = size(Nodes, 1);
                A = C\(K);
                
                
                ii = eye(3*nNodes, 3*nNodes);
                B = ii+dt*A;
                
                %[values] = eig(full(B), 'nobalance');
                [values] = eig(full(B));
                values = real(values);
                TheValue = max(abs(values));
                
                
                Kc = [Kc; k];
                Dtc = [Dtc; dt];
                Mc = [Mc; M];
                Hec = [Hec; he];
                Valc = [Valc; TheValue];
                
                %planeFitting?
                nn = length(Kc);
                A = [ones(nn,1) log10(Kc), log10(Mc), log10(Dtc), log10(Hec) ];
                b = log10(Valc);
                x = A\b
                hola = 1;
                toc
                
                relE = max(abs(A*x-b)./b)
            end
        end
    end
end