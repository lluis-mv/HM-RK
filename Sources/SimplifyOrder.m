function [NodOut, ElOut] = SimplifyOrder(Nodes,Elements)


nN = size(Elements,2);
El = Elements(:,1:nN/2);


nn = unique( El);

hola = 1;
ElOut = El;
for i = 1:length(nn)
    index = find( El == nn(i));
    ElOut(index) = i;
    NodOut(i,:) = Nodes(nn(i),:);
end
    
