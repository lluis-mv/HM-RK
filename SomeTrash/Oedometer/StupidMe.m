function [] = new_idea3K()

addpath('../')
% 1. Define the problem

k = 1E-8;
E = 8000;
h = 0.015;

ElementType = 'T3T3';

[A, he, M] = GetThisMatrix(h, k, E, 'T3T3');
[A2, ~, M] = GetThisMatrix(h, k, E, 'T6T3');
[A3, ~, M] = GetThisMatrix(h, k, E, 'T6T6');

ddt = 10.^linspace(-7,1);
figure(2)
i = 1;
for dt = ddt

    ii = eye(size(A));
    B = ii+dt*A;

    [values] = eig(full(B));
    values = real(values);
    TheValue(i) = max(abs(values));
    
    ii = eye(size(A2));
    B2 = ii+dt*A2;

    [values] = eig(full(B2));
    values = real(values);
    TheValue2(i) = max(abs(values));
    
    
    ii = eye(size(A3));
    B3 = ii+dt*A3;

    [values] = eig(full(B3));
    values = real(values);
    TheValue3(i) = max(abs(values));
    
    
    loglog(ddt(1:i), TheValue(1:i), 'b*-.', ddt(1:i), TheValue2(1:i), 'r*-.',  ddt(1:i), TheValue3(1:i), 'm*-.')
    xlabel('dt')
    ylabel('Spectral radii')
    yy = ylim();
    
    
    
    hold on
    xx = he^2/(6*M*k)
    plot([xx, xx], yy, 'r')
    xx = he^2/(200*M*k)
    plot([xx, xx], yy, 'b')
    xx = he^2/(2000*M*k)
    plot([xx, xx], yy, 'm')
    drawnow
    hold off

    i = i +1;
end


function [A, he, M] =  GetThisMatrix(he, k, E, ElementType)


RKMethod = 1;




dx = 0.1; dy = 0.1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';


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
figure(2312)
trimesh(Elements, Nodes(:,1), Nodes(:,2))
axis equal
axis off
drawnow
    
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



dt = 0;
[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, 0);


[C, K, X, f, fini] = ApplyBoundaryConditions(Nodes, Elements, C, K);
C = full(C);
K = full(K);
nNodes = size(Nodes, 1);
A = C\(K);

