function [] = ExampleOne()

close all
addpath('../')
% 1. Define the problem

T = 1E-7;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.075;

model = createpde(1);

dx = 0.4; dy = 1;

R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';


% First part. compute the eigenvalues
figure(1);
clf;
triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
drawnow
axis equal
axis off
print('ExampleOne-FemMesh', '-dpdf')

Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');


NStepsRef = 1;








for Stab = [1, 0]
    
    FirstTime = true;
    
    for j = 3:-1:1
        if ( j == 1)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
            Color = 'r';
        elseif (j == 2)
            ElementType = 'T3T3';
            Nodes = Nodes1;
            Elements = Elements1;
            Color = 'g';
        else
            ElementType = 'T6T3';
            Nodes = Nodes2;
            Elements = Elements2;
            Color = 'b';
        end
        
        
        dt = t/NStepsRef;
        
        if (j ~= 1)
             [U, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, NStepsRef, ElementType, Stab);
        else
             [U, GPInfo] = ComputeImplicitNonLinearProblemNodal(Nodes, Elements, CP, dt, NStepsRef, ElementType, Stab);
        end

        
        index = find(Nodes(:,1) == 0);
        WP = U(3*(index-1)+3);
        
        y = Nodes(index,2);
        [y, index2] = sort(y);
        WP = WP(index2);
        
        if ( all('T6T6' == ElementType))
            [WP, y] = CorrectInterpolation(WP, y);
        end
        
        figure(987+Stab)
        
        if ( FirstTime)
            FirstTime = false;
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, U);
            Xa = Xa(3*(index-1)+3);
            Xa = Xa(index2);
            plot(Xa, y, 'k', 'linewidth', 3, 'DisplayName', 'Reference');
            hold on
        end
        
        
        if ( all('T3T3' == ElementType))
            if ( j == 1)
                ElementType = ['NS-', ElementType];
            end
            plot(WP, y, [Color, '-.'] ,'DisplayName', ElementType, 'linewidth', 1.5)
            
        else
            plot(WP, y, Color, 'DisplayName', ElementType, 'linewidth', 1.5)
        end
        hold on
        ll = legend('location', 'best');
        set(ll, 'interpreter', 'latex')
        set(ll, 'location', 'best')
        xlabel('$p_w$ (kPa)', 'interpreter', 'latex')
        ylabel('$z$ (m)', 'interpreter', 'latex')
        set(gca, 'FontSize', 16)
        print(['ExampleOne-Solution-', num2str(Stab)], '-dpdf')
        
    end
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



function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% Other analytical solution...
nNodes = size(Nodes,1);
M = CP.M;
k = CP.k;
for nod = 1:nNodes
    xx = 1-Nodes(nod,2);
    TT = M * t*k;
    pw = 0;
    for m = 0:400
        aux = pi/2*(2*m+1);
        pw = pw + 2/aux * sin( aux * xx) * exp( - aux^2 * TT);
    end
    Xa(3*(nod-1)+3) = pw;
end


