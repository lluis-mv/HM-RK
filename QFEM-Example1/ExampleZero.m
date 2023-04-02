function [] = ExampleZero()

close all
addpath('../Sources')
% 1. Define the problem

T = 1E-7;


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
CP.k = 10000.0;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;





[Nodes, Elements] = CreateFEM(20);


% First part. compute the eigenvalues
figure(1);
clf;
ii = [1,2,3,4,1];
for elem = 1:size(Elements,1)
    Celem = Elements(elem,:);
    Xe = Nodes(Celem,:);
    plot(Xe(ii,1), Xe(ii,2), 'k')
    hold on
end

drawnow
axis equal
 axis off



Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size


dt = 1E-9
nSteps = 2;
       


[U, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, 'Q4Q4', 1);
index = find(Nodes(:,1) == 0);
WP = U(3*(index-1)+3);
y = Nodes(index,2);
[y, index2] = sort(y);
WP = WP(index2);

[U2, GPInfo] = ComputeImplicitNonLinearProblemNodalQuad(Nodes, Elements, CP, dt, nSteps, 'Q4Q4', 1);

WP2 = U2(3*(index-1)+3);
WP2 = WP2(index2);

[Xa] = ComputeAnalyticalSolution(Nodes, Elements, 'Q4Q4', dt*nSteps, CP, GPInfo, U);
Xa = Xa(3*(index-1)+3);
Xa = Xa(index2);

figure(4)
plot(Xa, y, 'k', 'linewidth', 3, 'DisplayName', 'Reference');
hold on
plot(WP, y, 'r', 'linewidth', 3)
plot(WP2, y, 'g', 'linewidth', 3)
return;

for Stab = [1, 0]
    
    FirstTime = true;
    
    for j = 3:-1:1
        if ( j == 1)
            ElementType = 'Q4Q4';
            Nodes = Nodes1;
            Elements = Elements1;
            Color = 'r';
        elseif (j == 2)
            ElementType = 'Q4Q4';
            Nodes = Nodes1;
            Elements = Elements1;
            Color = 'g';
        else
            ElementType = 'Q8Q4';
            Nodes = Nodes2;
            Elements = Elements2;
            Color = 'b';
        end
        
        
        dt = t/NStepsRef;
        
        if (j ~= 1)
             [U, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, NStepsRef, ElementType, Stab);
        else
             [U, GPInfo] = ComputeImplicitNonLinearProblemNodalQuad(Nodes, Elements, CP, dt, NStepsRef, ElementType, Stab);
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
        
        
        if ( all('Q4Q4' == ElementType))
            if ( j == 1)
                ElementType = ['NS-', ElementType];
            end
            if ( Stab > 0.5)
                ElementType = [ElementType, '. Stab'];
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
        print(['ExampleQne-Solution-', num2str(Stab)], '-dpdf')
        
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


function [Nodes, Elements] = CreateFEM( nx)

ny = floor( nx/0.4);

dx = 0.4/nx;
dy = 1/ny;

nx = nx+1;
ny = ny+1;

Elements = [];
for j = 1:ny-1
    for i = 1:nx-1
    
        Elements = [ Elements;
            (j-1)*nx+i, (j-1)*nx+i+1, (j)*nx+i+1, (j)*nx+i];
    end
end

Nodes = [];
for j = 1:ny
    for i = 1:nx
        Nodes = [ Nodes;
            (i-1)*dx, (j-1)*dy];
    end
end
