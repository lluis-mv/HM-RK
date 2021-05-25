function [] = ExampleOne()
figure(1); clf;
figure(26); clf;
addpath('../')
% 1. Define the problem

T = 1E-2;


CP.HydroMechanical = true;
CP.E = 1;
CP.nu = 0.0;
CP.k = 1;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.4;

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
he = mean(sqrt( mean([GPInfo(:,:).Weight])));

NSteps = 10.^linspace(0, 4, 10);
NSteps = 10.^linspace(0, 2, 10);
NSteps = floor(NSteps); NSteps = sort(NSteps);
NStepsRef = 1;


ddtt = t./NSteps;

ticks = [min(ddtt), max(ddtt)];
ticks = floor(log10(ticks));
ticks = 10.^(ticks(1):2:ticks(end));
ticks = [ticks, min(ddtt), max(ddtt)];
ticks = 10.^unique( log10(ticks));
ticks = unique(ticks);



for j = 2
    
    if ( j == 1)
        ElementType = 'T3T3';
        Nodes = Nodes1;
        Elements = Elements1;
        ThisNumber = 200;
    elseif (j == 2)
        ElementType = 'T6T3';
        Nodes = Nodes2;
        Elements = Elements2;
        ThisNumber = 6;
    else
        ElementType = 'T6T6';
        Nodes = Nodes2;
        Elements = Elements2;
        ThisNumber = 2000;
    end
    
    
    
    for RK = [1,3,8,-1, -2, -3, -10]%-1,1,3,8,-1, -2, -3]
%     for RK = [1,3,8,-11]%-1,1,3,8,-1, -2, -3]
        firstTime = true;
        i = 1;
        
        for nSteps = NSteps
            
            dt = t/nSteps;
            if ( RK > 0)
                Stab = 1;
                [U,GPInfo] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, Stab);
            elseif ( RK == -10)
                Stab = 0;
%                 [U, GPInfo] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType);
                [U, GPInfo] = ComputeImplicitLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, 123, Stab);

                Ui = U;
            elseif ( RK < 0)
                Stab = 0;
                [U,GPInfo] = ComputeRKILinearProblem3(Nodes, Elements, CP, dt, nSteps, ElementType, RK, Stab);
                Urk = U;
            end
            if ( firstTime)
                [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo,U);
            end
            [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes1, Elements1, GPInfo);
            
            
            
            i = i+1;
        end
        
        
        color = '';
        if ( RK == 1)
            color = 'r';
        elseif ( RK == 3)
            color = 'b';
        elseif (RK == 8)
            color = 'k';
        elseif ( RK == -1)
            color = 'c';
        elseif ( RK == -2)
            color = 'm';
        elseif ( RK == -5)
            color = 'g';
        elseif ( RK == -10)
            color = 'y';
        end
        
        figure(20+3*j)
        if (RK == 1)
            clf;
        end
        loglog( ddtt, L2, [color, '*-.'], 'DisplayName',  ['$L_2 p_w$ RK ', num2str(RK)]);
        hold on
        loglog( ddtt, L2U, [color, 'v-.'], 'DisplayName',  ['$L_2 u$ RK ', num2str(RK)]);
        hold on
        xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
        ylabel('Error norm', 'interpreter', 'latex');
        set(gca, 'FontSize', 14)
        drawnow
        
        
    end
    figure(20+3*j)
    
    drawnow
    yy = ylim();
    xx = (he)^2/(CP.k*CP.M*ThisNumber)*[1,1];
    plot(xx, yy, 'k-.', 'HandleVisibility','off')
    if ( yy(2) > 1E20)
        yy(2) = 1E20;
        yy(1) = 1E-10;
    end
    ylim(yy);
    drawnow
    xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
    xticks(ticks);
    legend('location', 'best', 'interpreter', 'latex', 'location', 'bestoutside')
    print(['ExampleOne-RK-', ElementType], '-dpdf')
    

end








function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 3, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


Xa = real(Xa);





function [A, nDirichlet, nNoDir] = GetAMatrix(Nodes, Elements, CP, dt, ElementType, RKMethod, AlphaStabM)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, ElementType);


[C, K ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, RKMethod, dt, false, AlphaStabM);


[C, K, X, fini, nDirichlet] = ApplyBoundaryConditions(Nodes, Elements, GPInfo, C, K);

A = C\(K);

nNoDir = [];
for i = 1:3*nNodes
    if ( any(i == nDirichlet))
        % do nothing
    else
        nNoDir = [nNoDir, i];
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
