function [] = OedometerMCC()

addpath('../')
addpath('../ModifiedCamClay/')


CP.E = 100;
CP.nu = 0.3;
CP.k = 1E-3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = 0.1;


eSize = 0.05;
eSize = 0.5;
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
% NSteps = [1,2];
i = 1;
errI = 0*NSteps;
errE = 0*NSteps;
errE2 = 0*NSteps;

err2E = 0*NSteps;
err2E2 = 0*NSteps;

for nSteps = NSteps
    dt = t/nSteps;
%      [U,GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps);
%      errI(i) = abs(GPInfo(end).StressNew(2)-11)
%     
    
    [U,GPInfo, errE2(i)] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 1E-8);
    errE(i) = abs(GPInfo(end).StressNew(2)-11)
    
   [U,GPInfo, err2E2(i)] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 2);
    err2E(i) = abs(GPInfo(end).StressNew(2)-11)
    
    
    figure(1)
    %loglog(NSteps, errI, 'r*-.', NSteps, errE, 'b*-.',  NSteps, errE2, 'g*-.')
    loglog( NSteps, errE, 'b*-.',  NSteps, errE2, 'r*-.', NSteps, err2E, 'c*-.',  NSteps, err2E2, 'm*-.')
    
    i = i+1;
end


