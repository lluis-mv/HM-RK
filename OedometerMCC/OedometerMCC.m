function [] = OedometerMCC()

addpath('../')
addpath('../ModifiedCamClay/')

t = 0.1;


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);




eSize = 0.015;
eSize = 0.015;
eSize = 0.2;
% eSize = 0.0075;
model = createpde(1);

dx = 1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';




NSteps = 10.^linspace(0, 3, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);

i = 1;
errI = 0*NSteps;
errE = 0*NSteps;
errE2 = 0*NSteps;

err2E = 0*NSteps;
err2E2 = 0*NSteps;

% NSteps = 200;
for nSteps = NSteps
    dt = t/nSteps;
%      [U,GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', 1);
%      errI(i) = abs(GPInfo(end).StressNew(2)+11)
    
    
    [U,GPInfo, errE2(i)] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', 1, false);
    errE(i) = abs(GPInfo(end).StressNew(2)+11)
    
    figure(1)
    loglog( NSteps, errE, 'b*-.',  NSteps, errE2, 'r*-.', NSteps, err2E, 'c*-.',  NSteps, err2E2, 'm*-.')
    drawnow
    
    [U,GPInfo, err2E2(i)] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', 1, true);
    err2E(i) = abs(GPInfo(end).StressNew(2)+11)
    
    
    figure(1)
    loglog( NSteps, errE, 'b*-.',  NSteps, errE2, 'r*-.', NSteps, err2E, 'c*-.',  NSteps, err2E2, 'm*-.')
    drawnow
    
    i = i+1;
end


