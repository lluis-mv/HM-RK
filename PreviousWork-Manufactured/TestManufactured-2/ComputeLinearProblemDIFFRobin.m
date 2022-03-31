% Solver for a linear problem

function [X, GPInfo, ThisInfo] = ComputeLinearProblemDIFFRobin(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, AlphaStabM)

if (nargout == 3)
    DoSomePostProcess = true;
else
    DoSomePostProcess = false;
end

if (nargin ==7)
    AlphaStabM = 1;
end



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);
[GPInfo] = InitializeConstitutiveLaw(CP, GPInfo);

[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);
[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

RealDofs = unique([unique( [GPInfo.dofsWP]), unique( [GPInfo.dofsU])]);
for i = 1:3*nNodes
    if ( any(i == RealDofs) == false)
        C(i,:) = 0;
        C(:,i) = 0;
        K(i,:) = 0;
        K(:,i) = 0;
        C(i,i) = 1;
        K(i,i) = 1;
    end
end

% First I have to find the nodes....
% Create a vector with the boundary conditions

nodesBC = find(Nodes(:,1) == min(Nodes(:,1)));

% now I have to create the structure or whatever
BC = [];
for el = 1:size(Elements,1)
    found1 = false;
    found2 = false;
    for i = 1:size(Elements,2)
        if ( any( Elements(el,i) == nodesBC) )
            if ( found1 == false)
                found1 = Elements(el,i);
            else
                found2 = Elements(el,i);
            end
        end
        if ( found1 ~= false && found2 ~= false)
            BC = [BC; found1, found2];
        end
    end
end


t = 1;
[f, uDirichlet] = ComputeForceVector(t, Nodes, Elements, GPInfo, CP);
f(nDirichlet) = 0;




Kbc = 0*K;
fbc = 0*f;
for nbc = 1:size(BC,1)
    indexWP = 3*(BC(nbc,:)-1)+3;
    h = norm( Nodes(BC(nbc,1),:)-Nodes(BC(nbc,2),:));
    Kbc(indexWP,indexWP) = Kbc(indexWP,indexWP) + h/6*[2,1;1,2];
    fbc(indexWP) = fbc(indexWP) + h*[1/2,1/2]';
end

K = K-Kbc;
f = 0*f-fbc;




K(nDirichlet, nDirichlet) = eye( length(nDirichlet), length(nDirichlet));
K = K(3:3:end,3:3:end);
f = f(3:3:end);
U = -K\f;
X(3:3:end) = U;
