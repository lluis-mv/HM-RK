[Nodes, Elements] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(Nodes, Elements);

figure(1); clf;
plot(NodesT(:,1), NodesT(:,2), 'k*')
axis equal
hold on

ind =  [1,4,2,5,3,6,1];
for el = 1:size(ElementsT,1)
    Cel = ElementsT(el,:);
    plot( NodesT( Cel(ind), 1 ), NodesT( Cel(ind), 2));
end