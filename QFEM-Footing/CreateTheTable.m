function [] = CreateTheTable()
addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-12;
CP.k = 1E-12;
CP.Elastic = false;
CP.MCC = 2;




ESIZE = [2:11];


i = 1;




for eSize = ESIZE

    


   [Nodes, Elements] = ReadTheMesh(['Mesh', num2str(eSize), '.msh']);
   [Nodes1, Elements1] = SimplifyOrder(Nodes,Elements);



    [GPInfo] = ComputeElementalMatrices(Nodes1, Elements1, CP, 'Q4Q4');
    he = mean(sqrt([GPInfo(:,:).Weight]));

    eSizeAxis(i) = he;
    nElemens(i) = size(Elements1,1);
    nPatch(i) = size(Nodes1,1);
    nDofs(i) = size(Nodes1,1)*3;
    nDofsQ(i) = size(Nodes,1)*3;
    i = i+1;

    line = [ num2str(he), ' & ', num2str( size(Elements1,1)) ,' & ', num2str( size(Nodes1,1)), ' & ', ...
        num2str( size(Nodes1,1)*3), ' & ' num2str( size(Nodes,1)*3)]
end




hola = 1;