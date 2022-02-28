function [NodesT, ElementsT] = ConvertToTriangles(Nodes, Elements)

nNodEnd = size(Nodes, 1);
nElem = size(Elements,1);

NewNodes = [];
tAux1 = [1,2,4,5];
tAux2 = [2,3,4,6,7];
ElementsT = [];
for el = 1:nElem
    Cel = Elements(el,:);
    
    nNodEnd = nNodEnd+1;
    ElementsT = [ElementsT; 
        Cel(tAux1), nNodEnd, Cel(8);
        Cel(tAux2), nNodEnd];
    NewNodes = [NewNodes;
        0.5*( Nodes(Cel(2),:)+Nodes(Cel(4),:))];
    
    
end
NodesT = [Nodes; NewNodes];