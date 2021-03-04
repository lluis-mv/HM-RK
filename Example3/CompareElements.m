function [] = CompareThis()

figure(212); hold off; clf;
addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-6;
CP.Elastic = false;
CP.MCC = true;

eSize = [0.45];




LinearElastic = false;
if ( LinearElastic)
    CP.Elastic = true;
    CP.MCC = false;
end


% % if (Elem == 1)
% %     ElementType = 'T3T3';
% % elseif (Elem == 2)
% %     ElementType = 'T6T3';
% % else
% %     ElementType = 'T6T6';
% % end

model = createpde(1);


R1 = [3,5, 0, 1, 4, 4, 0, 0, 0, 0, -4, -4]';

% R1 = [3,4, 0, 0.4, 0.4, 0, 0, 0, -1, -1]';

g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize);
Nodes = mesh.Nodes';
Elements = mesh.Elements';



% Estimate the element size
modelA = createpde(1);
geometryFromEdges(modelA, g);
mesha = generateMesh(modelA, 'Hmax', eSize, 'GeometricOrder','linear');
Nodesa = mesha.Nodes';
Elementsa = mesha.Elements';
[GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));


RK = 1;

Perme = 10.^(-2);
CP.k = Perme;

nSteps = 50;
dt = 0.15/nSteps;

if ( LinearElastic)
    [U, GPInfo,  information] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', RK, 1);
else
    [U, GPInfo, rrr,  information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3', RK, 1, false);
end
figure(555)
pdeplot(model,'XYData',U(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

FF = [information.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), 'DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on



if ( LinearElastic)
    [U2, GPInfo2,  information2] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T6', RK, 1);
else
     [U2, GPInfo2, rrr2,  information2] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, 'T6T6', RK, 1, false);
end
figure(556)
pdeplot(model,'XYData',U2(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;



FF = [information2.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), 'DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on



model1 = createpde(1);




g = decsg(R1);
geometryFromEdges(model1, g);


mesh1 = generateMesh(model1, 'Hmax', eSize, 'GeometricOrder','linear');



Nodes1 = mesh1.Nodes';
Elements1 = mesh1.Elements';



if ( LinearElastic)
    [U1, GPInfo1,  information1] = ComputeLinearProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', RK, 1);
else
    [U1, GPInfo1, rrr1,  information1] = ComputeNLProblem(Nodes1, Elements1, CP, dt, nSteps, 'T3T3', RK, 1, false);
end
figure(559)
pdeplot(model1,'XYData',U1(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;



FF = [information2.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), '-.', 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), '-.','DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on





[Ur, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T6T3');
figure(88)
pdeplot(model,'XYData',Ur(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

figure(100)
pdeplot(model,'XYData',Ur(3:3:end)-U(3:3:end),'ColorMap','jet', 'Mesh','on');
drawnow;

FF = [information.F];
figure(212)
subplot(2,1,1)
plot( [information.t], FF(1:2:end), 'DisplayName', ['k=', num2str(Perme)])
hold on
subplot(2,1,2)
plot( [information.t], FF(2:2:end), 'DisplayName', ['k=', num2str(Perme)])
legend('location', 'best')
hold on

[A, nDir, nnoDir] = GetAMatrix(Nodes1, Elements1, CP, dt, 'T3T3', 1, 1);
nNodes = size(Nodes1, 1);
ii = eye(3*nNodes, 3*nNodes);

B = ii + dt*A;
B = B(nnoDir, nnoDir);

values = eig(full(B), 'nobalance');
values = abs(values);
minval= min(values);
maxval = max(values);
hola = 1;


    
maxx = 0;
for i = [88, 555, 559, 556]
    figure(i)
    xx = caxis();
    maxx = maxx+xx(2);
end

for i = [88, 555, 559, 556]
    figure(i)
    caxis([0, maxx/4])
end
    
            



function [A, nDirichlet, nNoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, RKMethod, AlphaStabM)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);


[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);


[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

A = C\(K);

nNoDir = [];
for i = 1:3*nNodes
    if ( any(i == nDirichlet))
        % do nothing
    else
        nNoDir = [nNoDir, i];
    end
end

