close all
addpath('../')

[Nodes, Elements] = ReadTheMesh('ThisMesh.msh');
[Nodes, Elements] = ConvertToTriangles(Nodes, Elements);
figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-1'], '-dpdf')


[Nodes, Elements] = ReadTheMesh('ThisMeshT.msh');
figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-2'], '-dpdf')

[Nodes, Elements] = ReadTheMesh('MeshUglyT.msh');
figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-3'], '-dpdf')



[Nodes, Elements] = ReadTheMesh('ThisMesh.msh');

figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-4'], '-dpdf')


[Nodes, Elements] = ReadTheMesh('ThisMeshQ.msh');
figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-5'], '-dpdf')

[Nodes, Elements] = ReadTheMesh('MeshUglyQ.msh');
figure(380); clf;
PlotMesh(Nodes, Elements);
axis off
axis equal
print(['ExampleOne-FemMesh-6'], '-dpdf')

