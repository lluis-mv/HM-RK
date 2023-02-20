% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, AlphaStabM)


if (nargout == 3)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 1, dt, true, AlphaStabM);

[C, K, X, fini, nDir] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);


if ( ElementType(1) == 'T')
PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, 0, true, ['ThisIProblem-', ElementType]);
end
if ( DoSomePostProcess )
    ThisInfo = DoThisPostProcess( 0, Nodes, Elements, GPInfo, X, CP);
end

A = C\(K);
ii = eye(3*nNodes, 3*nNodes);
B = ii-dt*A;
B = inv(B);


B(nDir,nDir) = eye(length(nDir));
invC = inv(C);
invCfini = (C\fini);
invC(nDir,nDir) = 0;
invCfini(nDir) = 0;

[f] = ComputeForceVector(dt, Nodes, Elements, GPInfo);
X = B*(X + invCfini + (1/nSteps)*invC*f);

for i = 2:nSteps
    [f] = ComputeForceVector(i*dt, Nodes, Elements, GPInfo);
    
    X = B*(X +  (1/nSteps)*invC*f);
	if ( DoSomePostProcess )
        GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false, dt);
        GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);
        ThisInfo = DoThisPostProcess( i*dt, Nodes, Elements, GPInfo, X, CP, ThisInfo);
    end
end


GPInfo = EvaluateConstitutiveLaw(CP, GPInfo, X, Elements, false, dt);
GPInfo = FinalizeConstitutiveLaw(CP, GPInfo);

if ( ElementType(1) == 'T')
PostProcessResults(CP.HydroMechanical, Nodes, Elements, X, GPInfo, dt*nSteps, false, ['ThisIProblem-', ElementType]);
end
