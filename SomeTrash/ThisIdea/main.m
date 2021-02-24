function [] = ExampleOneBis()

close all

% 1. Define the problem



ESIZE = [0.2, 0.1, 0.075, 0.05, 0.025, 0.02];
ESIZE = [0.2, 0.18, 0.15, 0.1, 0.075, 0.06, 0.05, 0.04, 0.03, 0.025, 0.02];
figure(50); clf;

for Elem = 1:2
    
    esizeAxis = ESIZE;
    i = 1;
    for eSize = ESIZE
        
        if (Elem == 1)
            ElementType = 'T3';
            ThisNumber = 200;
        elseif (Elem == 2)
            ElementType = 'T6';
            ThisNumber = 6;
        end
        figure(20)
        dx = 1; dy = 1;
        model = createpde(1);
        
        R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
        g = decsg(R1);
        geometryFromEdges(model, g);
        pdegplot(model,'EdgeLabels','on'); 
        axis equal
        
        
        applyBoundaryCondition(model,'dirichlet','Edge',2,'u',1);
        applyBoundaryCondition(model,'dirichlet','Edge',4,'u',0);
        fun =  @(x,y) -12*(x.x).^2;
        specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',fun);
        
        if ( Elem == 1)
            mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        else
            mesh = generateMesh(model, 'Hmax', eSize);
        end
        
        Nodes = mesh.Nodes';
        Elements = mesh.Elements';
        
        results = solvepde(model);
        u = results.NodalSolution;
        pdeplot(model,'XYData',u)
        title('Numerical Solution');
        xlabel('x')
        ylabel('y')
        
        
        
        
        esizeAxis(i)=eSize;
        
        
        
        
        [L2(i), LInf(i)] = ComputeErrorNorms(u, Nodes, Elements);
        
        
        
        figure(30)
        
        loglog( esizeAxis(1:i), L2(1:i), 'k*-.',  esizeAxis(1:i), LInf(1:i), 'g*-.')
        hold on
        
        i = i+1;
        
        
    end
    
    % Slope, que jode....
    clear Slope1;
    clear Slope2
    for i = 2:length(esizeAxis)
        Slope1(i) = log10(LInf(i)/LInf(i-1)) / log10(esizeAxis(i)/esizeAxis(i-1));
        Slope2(i) = log10(L2(i)/L2(i-1)) / log10(esizeAxis(i)/esizeAxis(i-1));
    end
    Slope1
    Slope2
    
    figure(50)
    
    loglog( esizeAxis, LInf, 'g*-.',  esizeAxis, L2, 'k*-.')
    hold on
    xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
    ylabel('Error norm', 'interpreter', 'latex');
    set(gca, 'FontSize', 14)
    drawnow
end









function [L2, LInf] = ComputeErrorNorms(X, Nodes, Elements)

[Xa] = ComputeAnalyticalSolution(Nodes, Elements, X);

nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


for i = 1:nNodes
    LInf= max(abs(X-Xa));
end



alfa = 2/3; beta = 1/6;
N1 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 1/6;
N2 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 2/3;
N3 = [ 1 - alfa - beta, alfa,  beta];

wa = 0.054975871827661;
wb = 0.1116907948390055;
Na1 = 0.816847572980459;
Nb1 = 0.108103018168070;
Na2 = 0.091576213509771;
Nb2 = 0.445948490915965;



auxK = [Na2, Na2, wa;
    Na1, Na2, wa;
    Na2, Na1, wa;
    Nb2, Nb2, wb ;
    Nb1, Nb2, wb ;
    Nb2, Nb1, wb ];

al = auxK(:,1)';
be = auxK(:,2)';
w = auxK(:,3)'/sum(auxK(:,3)');


L2 = 0;

[GPInfo] = ComputeElementalMatrices(Nodes, Elements);
for el = 1:nElements
    Cel = Elements(el,:);
    dofs = Cel;
    
    
    wA = sum([GPInfo(el,:).Weight]);
    
    for gp = 1:length(w)
        [Np] = GetShapeFunctions( al(gp), be(gp), length(dofs));
        L2 = L2 + wA*w(gp)*abs( Np * ( X(dofs)-Xa(dofs)));
    end
    
end

L2 = L2/sum([GPInfo.Weight]);



if (L2 < 1E-15)
    L2 = rand*1E-15;
end

if (LInf < 1E-15)
    LInf = rand*1E-15;
end



function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, Xnum)
Xa = 0*Xnum;
nNodes = size(Nodes, 1);
nElements = size(Elements, 1);

for nod = 1:nNodes
    Xa(nod) = Nodes(nod,1).^4;
end

function [Np] = GetShapeFunctions( alfa, beta,nP)

ndim = 2;
nU = 6;
if (nU == 6)
    Nsmall =  [ 1 - alfa - beta; alfa;  beta];
    
    Nu = (zeros(ndim, 3*ndim));
    for i = 1:3
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
elseif (nU == 12)
    Nsmall =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)];
    
    
    Nu = (zeros(ndim, 6*ndim));
    for i = 1:6
        for dd = 1:2
            Nu(dd, ndim*(i-1)+dd) = Nsmall(i);
        end
    end
end


if ( nP == 3)
    Np =  [ 1 - alfa - beta; alfa;  beta]';
elseif ( nP == 6)
    Np =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
        alfa*(2*alfa-1);
        beta*(2*beta-1);
        4*alfa*(1-alfa-beta);
        4*alfa*beta;
        4*beta*(1-alfa-beta)]';
end






function [p, y] = CorrectInterpolation(p1, y1)

alfa = linspace(0, 1);
y = [];
p = [];

for ind = 1:2:length(y1)-1
    pp = p1(ind:ind+2);
    yy = y1(ind:ind+2);
    
    N = [(1 - alfa).*(1-2*alfa);
        4*(1-alfa).*alfa;
        alfa.*(2*alfa-1)];
    y = [y; N'*yy];
    p = [p; N'*pp];
end