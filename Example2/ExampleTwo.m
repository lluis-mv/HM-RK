function [] = ExampleOne()

figure(31); hold off
figure(32); hold off
figure(33); hold off

addpath('../')
% 1. Define the problem

T = 1;


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;

nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = 1;


eSize = 0.1;

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
% print('ExampleOne-FemMesh', '-dpdf')

Nodes1 = Nodes;
Elements1 = Elements;

% Estimate the element size
[GPInfo] = ComputeElementalMatrices(Nodes, Elements, CP, 'T3T3');
he = mean(sqrt( mean([GPInfo(:,:).Weight])));

NSteps = 10.^linspace(0, 2, 5);
NSteps = floor(NSteps); NSteps = sort(NSteps);
NStepsRef = 1;


ddtt = t./NSteps;

ticks = [min(ddtt), max(ddtt)];
ticks = floor(log10(ticks));
ticks = 10.^(ticks(1):2:ticks(end));
ticks = [ticks, min(ddtt), max(ddtt)];
ticks = 10.^unique( log10(ticks));
ticks = unique(ticks);




for j = 1:3
    % Now the same to compute norms...
    for drift = [false, true]
        for RKMethod = [1:7]
            
            
            
            if ( j == 1)
                ElementType = 'T3T3';
                Nodes = Nodes1;
                Elements = Elements1;
                ThisNumber = 200;
                CP.k = 1E-8;
            elseif (j == 2)
                ElementType = 'T6T3';
                Nodes = Nodes2;
                Elements = Elements2;
                ThisNumber = 6;
                CP.k = 1E-7;
            else
                ElementType = 'T6T6';
                Nodes = Nodes2;
                Elements = Elements2;
                ThisNumber = 2000;
                CP.k = 1E+3;
            end
            
            
            ElementType = 'T6T6';
            Nodes = Nodes2;
            Elements = Elements2;
            ThisNumber = 6;
            CP.k = 1E-4;
            
            
            i = 1;
            for nSteps = NSteps
                
                dt = 1/nSteps;
                
                [U, GPInfo, RES] = ComputeThisNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, 1, drift);
                [U, GPInfo, RES] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RKMethod, 1, drift);
                %[U2, GPInfo2] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType);  
                
                errResidu(i) = RES;
                
                figure(30)
                loglog( ddtt(1:i), errResidu(1:i), 'r*-.')
                hold on
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex');
                set(gca, 'FontSize', 14)
                drawnow
                
                ll = legend('$L_2 p_w$', '$L_2 u$', '$L_\infty p_w$', '$L_\infty u$', 'location', 'best');
                set(ll, 'interpreter', 'latex')
                drawnow
                hold off
                
                i = i+1;
            end
            
            figure(30+j)
            merda = loglog( ddtt, errResidu, '*-.');
            if ( drift)
                set(merda, 'DisplayName', ['RK-', num2str(RKMethod), ' drift']);
            else
                set(merda, 'DisplayName', ['RK-', num2str(RKMethod)]);
            end
            hold on
            xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
            ylabel('Error norm', 'interpreter', 'latex');
            set(gca, 'FontSize', 14)
            drawnow
            ll = legend('location', 'best');
            set(ll, 'interpreter', 'latex')
            drawnow
            xlim([0.9999*min(ddtt), 1.0001*max(ddtt)])
            xticks(ticks);
            print(['ExampleTwo-', num2str(j)], '-dpdf')
            
            
            %         hold off
        end
    end
end




