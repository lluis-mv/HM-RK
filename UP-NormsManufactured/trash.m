eSize = 0.1;
dx = 0.3; dy = 1;
model = createpde(1);





R1 = [3,5, 0, 1, 2, 2, 0, 0, 0.00001, 0, -2, -2]';
    

g = decsg(R1);
geometryFromEdges(model, g);


mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');


Nodes = mesh.Nodes';
Elements = mesh.Elements';

Stab = 1;
% First part. compute the eigenvalues
figure(1);
clf;
triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
drawnow
axis equal
axis off

