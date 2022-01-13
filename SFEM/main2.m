% an error something has to be programed for the splitting. Hmax should be
% the previous step Hmax and not the previous iteration one. It requires an
% adaptive substepping in time and reject steps
% Something to do displacement load curves has to be done to see the
% softening and whatever. xD. besis

function []=main()
figure(1)
% mesh
model = createpde(1);
R1 = [3,4,-3,3,3,-3,1,1,-1,-1]';
dx = 1; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]'; 
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', 0.15, 'GeometricOrder','linear');
pdeplot(model)

X = mesh.Nodes';
C = mesh.Elements';
nNodes = size(X,1);
nElem = size(C, 1);
nDim = 2;

U = zeros( nNodes*nDim, 1);

%pause()
% dirichlet boundary conditions
nodesBottom = find( X(:, 2) == min(X(:,2)) )';
nodesTop = find( X(:, 2) == max(X(:, 2)) )';
nodesLeft = find( X(:, 1) == min(X(:, 1)))';
nodesRight = find( X(:, 1) == max(X(:, 1)))';

%nDirichlet = [ 2.*(nodesBottom)-1  2.*(nodesBottom ) 2.*(nodesTop)-1  2.*(nodesTop),    2.*(nodesLeft)-1, 2.*(nodesRight)-1 ];
%UDirichlet = [ zeros(1, 2*length(nodesBottom))  0.0*ones(1, length(nodesTop)) 0.01*ones(1,length(nodesTop)), zeros(1, length(nodesLeft)), zeros(1, length(nodesRight)) ];

BC.nDirichlet = [ 2.*(nodesBottom)-1,  2.*(nodesBottom ),    2.*(nodesLeft)-1, 2.*(nodesRight)-1 ];
BC.UDirichlet = [ zeros(1, 2*length(nodesBottom)), zeros(1, length(nodesLeft)), zeros(1, length(nodesRight)) ];

BC.nLoad = nodesTop;
BC.nLoad = CreateLoadSegments( X, C, nodesTop);
BC.Load = [0, -1];


ElementalMatrices = ConstructElementalBMatrix( X, C);
NodalMatrices = ConstructNodalBMatrix( X, C, ElementalMatrices);

iU1 = [1:2:2*nNodes-1];
iU2 = [2:2:2*nNodes];



for i = 1:1
    [UNew, STRESS] = ComputeMechanicalProblem(X, C, U, BC, ElementalMatrices, NodalMatrices);
    U = UNew;
    figure(1)
    pdeplot(model, 'XYData', STRESS(:,1))
    figure(2)
    pdeplot(model, 'XYData', STRESS(:,2))
    figure(3)
    pdeplot(model, 'XYData', STRESS(:,3))
    figure(11)
    pdeplot(model, 'XYData', UNew(iU1))
    figure(12)
    pdeplot(model, 'XYData', UNew(iU2))
    figure(21)
    triplot( C, X(:, 1), X(:, 2))
    hold on
    quiver( X(:,1),X(:,2), UNew(iU1), UNew(iU2) );
    axis equal;
    hold off;
end




function [UNew, STRESS] = ComputeMechanicalProblem( X, C, U, BC, EM, NM)

%putDirichlet conditions
U(BC.nDirichlet) = BC.UDirichlet;

%Construct the matrices
[K, f] = EnsambleMechanicalProblem( X, C, U, BC, EM, NM);

dUNew = K\f;
UNew = U - dUNew;
%putDirichlet conditions (again)
UNew(BC.nDirichlet) = BC.UDirichlet;

% compute stresses and so on
% UNew = 0*UNew;
% UNew([2:2:end]) = X(:,2);
[K, f, STRESS] = EnsambleMechanicalProblem( X, C, UNew, BC, EM, NM);


function [K, f, STRESS] = EnsambleMechanicalProblem(X, C, U, BC, EM, NM)

nNodes = size(X,1);
nDim = 2;
nElem = size(C,1);

K = sparse( 2*nNodes, 2*nNodes);
f = zeros(2*nNodes, 1);
G = 800;
Poisson = 0.0;


Bulk = 2*G/3;
Bulk = 2 * G * ( 1 + Poisson) / 3 / ( 1- 2 * Poisson);
STRESS = zeros( size(NM,2), 3);
E = 3 * Bulk * ( 1- 2 * Poisson);


D = E / (1+Poisson) / (1-2*Poisson)  * [(1-Poisson), Poisson, 0; Poisson, (1-Poisson), 0; 0 0 (1-2*Poisson)/2];


for n = 1:size(NM,2)
    Ce = NM(n).NeigNodes;
    Xe = X(Ce,:);
    index = zeros(2*length(Ce),1);
    
    for i = 1:length(Ce)
        index(2*(i-1) +1 ) = 2*(Ce(i)-1) + 1;
        index(2*(i-1) +2 ) = 2*(Ce(i)-1) + 2;
    end
    ii = index;
    
    UPatch = U(index);
    BPatch = NM(n).B;
    Strain = BPatch' * UPatch;
    Stress = D * Strain;
    
    weight = NM(n).Area;
    RHS = BPatch * Stress * weight;
    LHS = BPatch * D * BPatch' * weight;
    
    %[RHS2, LHS2, stress2, psi2 ] = mechanicalProblem( G,  Bulk, 0, 0,0, U(ii(1)), U(ii(2)), U(ii(3)), U(ii(4)), U(ii(5)), U(ii(6)), ...
    %    Xe(1,1), Xe(1,2), Xe(2,1), Xe(2,2), Xe(3,1), Xe(3,2));
    %RHS = RHS2; LHS = LHS2; stress = stress2;
    K(index,index) = K(index,index)+LHS;
    f(index) = f(index) + RHS;
    
    STRESS(n,:) = Stress;
    
    
end

% Ensamble loads
for n = 1:size(BC.nLoad,1)
    Ce = BC.nLoad(n, :);
    Xe = X(Ce,:);
    
    dx = norm( Xe(1,:)-Xe(2,:))/2;
    
    index = zeros(2*length(Ce),1);
    for i = 1:length(Ce)
        index(2*(i-1) +1 ) = 2*(Ce(i)-1) + 1;
        index(2*(i-1) +2 ) = 2*(Ce(i)-1) + 2;
    end
    
    f(index) = f(index) - dx * [BC.Load, BC.Load]';
end

% put dirichlet boundary conditions
K(BC.nDirichlet,:) = 0;
f(BC.nDirichlet) = 0;
K(BC.nDirichlet, BC.nDirichlet) =  eye( size(BC.nDirichlet,2));
disp( ['NORM: ' num2str( norm(f)) ] )



function [ElementalMatrices] = ConstructElementalBMatrix( X, C)

% Construct and save the B matrix of each element
Nsmall_chi = [-1 -1; 1 0; 0 1];
nNodes = 3;
for el = 1:size(C,1)
    Celem = C(el,:);
    Xelem = X(Celem,:);
    J = (zeros(2,2));
    for j = 1:nNodes
        J = J + [ Nsmall_chi(j,1) * Xelem( j,1), Nsmall_chi(j,1) * Xelem( j,2);
            Nsmall_chi(j,2) * Xelem( j,1), Nsmall_chi(j,2) * Xelem( j,2)];
    end
    dN_dX = J \ Nsmall_chi';
    dN_dX = dN_dX';
    B = zeros(6,3);
    for j = 1:nNodes
        B(2*(j-1)+1, 1) = dN_dX(j,1);
        B(2*(j-1)+2, 2) = dN_dX(j,2);
        B(2*(j-1)+1, 3) = dN_dX(j,2);
        B(2*(j-1)+2, 3) = dN_dX(j,1);
    end
    ElementalMatrices(el).B = B;
    
    Area = [1 1 1;
    Xelem(1,1) Xelem(2,1) Xelem(3,1);
    Xelem(1,2) Xelem(2,2) Xelem(3,2)];

    ElementalMatrices(el).Area = det(Area)/2;
end




function NM = ConstructNodalBMatrix( X, C, BMatrices)

nNodes = size(X,1);

for i = 1:nNodes
    [candidates, trash] = find( C == i);
    NM(i).NeigElem = candidates;
    NM(i).NeigNodes = unique(C(candidates,:));
    NM(i).Area = sum( [BMatrices(NM(i).NeigElem).Area])/3;
end

for i = 1:nNodes
    Bnod = zeros( 2*length( NM(i).NeigNodes), 3);
    areas = 0;
    for neigElem = NM(i).NeigElem'
        B = BMatrices(neigElem).B * (BMatrices(neigElem).Area / 3) / NM(i).Area;
        Celem = C(neigElem,:);
        ind = [];
        for j = 1:3
            abc = find( NM(i).NeigNodes == Celem(j));
            ind = [ind, 2*(abc)-1, 2*abc];
        end
        Bnod(ind,:) = Bnod(ind,:) + B;    
    end
    
    NM(i).B = Bnod;
end


return;
clear NM
for i = 1:size(C,1)
    NM(i).NeigNodes = C(i,:);
    NM(i).B= BMatrices(i).B;
    NM(i).Area = BMatrices(i).Area;
end


function [segments] = CreateLoadSegments( X, C, nodesList)
segments = [];
nNodes = length(nodesList);
nodesList = sort(nodesList);
for i = 1:nNodes
    node = nodesList(i)
    [candidates, trash1] = find( C == node);
    CC = C(candidates,:);
    for j = i+1:nNodes
        otherNode = nodesList(j);
        [candidates, trash2] = find( CC == otherNode);
        if ( length(candidates) > 0)
            segments = [segments; node, otherNode];
        end
    end
    
    
end    
    