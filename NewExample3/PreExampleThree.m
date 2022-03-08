function [] = PreExampleThree()



addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = false;
CP.MCC = true;

ESIZE = 0.35;



RKReference = 8;
RKMethods = [8, 2, 12, 4,14, 6, 16];

for Elem = [1,2,3]
    
    esizeAxis = ESIZE;
    i = 1;
    for eSize = ESIZE
        
        if (Elem == 1)
            ElementType = 'T3T3';
        elseif (Elem == 2)
            ElementType = 'T6T3';
        else
            ElementType = 'T6T6';
        end
        
        model = createpde(1);
        
        
        R1 = [3,5, 0, 1, 3, 3, 0, 0, 0, 0, -3, -3]';
        
        g = decsg(R1);
        geometryFromEdges(model, g);
        
        if ( Elem == 1)
            mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        else
            mesh = generateMesh(model, 'Hmax', eSize);
        end
        
        Nodes = mesh.Nodes';
        Elements = mesh.Elements';
        
        figure(1);
        clf;
        if ( Elem == 1)
            triplot(Elements, Nodes(:,1), Nodes(:,2), 'k');
        end
        drawnow
        axis equal
        axis off
        
        
        % Estimate the element size
        
        mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodesa = mesha.Nodes';
        Elementsa = mesha.Elements';
        [GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
        he = mean(sqrt( mean([GPInfo(:,:).Weight])));
        
        Nadim = 20;
        
        NSteps = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8];
        for j = 1:length(NSteps)
            for RK = RKMethods
                nSteps = NSteps(j);
                dt = 0.15/nSteps;
                
                if ( RK > 10)
                    RK = RK-10;
                    tic;
                    [U,GPInfo, rrr, information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1, false);
                    disp('caseImplex')
                    toc
                    RK = RK+10;
                    timing = toc;
                else
                    tic;
                    [U,GPInfo, rrr, information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1, true);
                    disp('caseNoImplex')
                    toc
                    timing = toc;
                end
                
                
                ThisInfo(RK,j).t = [information.t];
                ThisInfo(RK,j).F = [information.F];
                ThisInfo(RK,j).F = ThisInfo(RK,j).F(1:2:end);
                
                if ( RK == RKReference)
                    Nadim = ThisInfo(RK,j).F(end);
                end
                
                RES(RK,j) = rrr;
                ddtt(RK,j) = dt;
                N(RK,j) = ThisInfo(RK,j).F(end);
                Time(RK,j) = timing;
                if ( size(N,1) >= 8)
                    N(8,:) = 0;
                end
                ind = find(N == 0);
                N(ind) = nan; ddtt(ind) = nan; Time(ind) = nan; RES(ind) = nan;
                
                if ( j < 2 )
                    continue;
                end
                
                figure(2105); clf
                for jj = 1:size(RES,1)
                    if ( all(isnan(RES(jj,:))) )
                        continue;
                    end
                    merda = SelectColor(jj);
                    loglog(ddtt(jj,:), RES(jj,:), [merda])
                    hold on
                    
                end
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Residual', 'interpreter', 'latex')
                
                grid minor
                drawnow
                hold off
                
                figure(2106); clf
                for jj = 1:size(RES,1)
                    if ( all(isnan(N(jj,:))) )
                        continue;
                    end
                    merda = SelectColor(jj);
                    loglog(ddtt(jj,:), N(jj,:), [merda])
                    hold on
                end
                grid minor
                drawnow
                hold off
                
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Footing load', 'interpreter', 'latex')
                
                figure(2107); clf
                
                for jj = 1:size(RES,1)
                    if ( all(isnan(RES(jj,:))) )
                        continue;
                    end
                    merda = SelectColor(jj);
                    plotThis = abs(N(jj,:)-Nadim);
                    if ( jj ~= RKReference)
                        index = find( plotThis == 0);
                        plotThis(index) = 1E-14*(1+rand(size(index)));
                    end
                    loglog(ddtt(jj,:), plotThis, [merda])
                    hold on
                end
                
                
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                
                figure(2108); clf
                ind = find(Time == 0);
                Time(ind) = nan;
                for jj = 1:size(RES,1)
                    if ( all(isnan(RES(jj,:))) )
                        continue;
                    end
                    merda = SelectColor(jj);
                    plotThis = abs(N(jj,:)-Nadim);
                    if ( jj ~= RKReference)
                        index = find( plotThis == 0);
                        plotThis(index) = 1E-14*(1+rand(size(index)));
                    end
                    plotThis(index) = 1E-14*(1+rand(size(index)));
                    loglog(Time(jj,:), plotThis, [merda])
                    hold on
                end
                
                
                xlabel('Computational cost (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                for ff = [2105:1:2108]
                    figure(ff);
                    ll = legend('2','4','6','2-E','4-E','6-E','location','bestoutside');
                    set(ll, 'interpreter', 'latex')
                end
            end
            % Printing....
            figure(2105)
            pause(1)
            print(['ExampleThree-Prev-Residual-', num2str(Elem)], '-dpdf')
            figure(2106)
            pause(1)
            print(['ExampleThree-Prev-Bearing-', num2str(Elem)], '-dpdf')
            figure(2107)
            pause(1)
            print(['ExampleThree-Prev-Error-', num2str(Elem)], '-dpdf')
            figure(2108)
            pause(1)
            print(['ExampleThree-Prev-TimeError-', num2str(Elem)], '-dpdf')
        end
        
        clear Time;
        clear RES;
        clear ddtt;
        clear N;
        clear Time
        
    end
    
end


function [merda ] = SelectColor(j)

switch j
    case 2
        merda = 'r-^';
    case 4
        merda = 'g-o';
    case 6
        merda = 'b-s';
    case 12
        merda = 'r-.^';
    case 14
        merda = 'g-.o';
    case 16
        merda = 'b-.s';
end
