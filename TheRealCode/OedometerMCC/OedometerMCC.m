function [] = OedometerMCC()

addpath('../')
% 1. Define the problem




CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-2;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = 0.1;


eSize = 0.05;


model = createpde(1);

dx = 0.1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

% figure(1)
% pdeplot(model)
% drawnow


Nodes = mesh.Nodes';
Elements = mesh.Elements';





NSteps = 10.^linspace(0, 3, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);

i = 1;
err = 0*NSteps;
for nSteps = NSteps
    dt = t/nSteps;
    [U,GPInfo] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps);
    err(i) = abs(GPInfo(end).StressNew(2)-30)
    loglog(NSteps, err, '*-.')
    i = i+1;
end


