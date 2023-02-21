function [] = PreExampleThree()



addpath('../Sources')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = false;
CP.MCC = true;

% Modified Cam clay parameters
CP.kappa = 0.01;
CP.lambda = 0.1;
CP.M_MCC = 1;
CP.nu = 0.3;
CP.u = 20;

ESIZE = 0.35;



RKReference = 8;
RKMethods = [8, 2, 12, 4,14, 6, 16];

[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);

for Elem = [2,1]
    
    esizeAxis = ESIZE;
    i = 1;
    for eSize = ESIZE
        
        
        if (Elem == 1)
            ElementType = 'T6T3';
            Elements = ElementsT;
            Nodes = NodesT;
        elseif (Elem == 2)
            ElementType = 'Q8Q4';
            Elements = ElementsQ;
            Nodes = NodesQ;
        end
        
        
        
                
        
        % Estimate the element size
        
        
        
        Nadim = 20;
        
        NSteps = [10, 2^4, 2^5, 2^6, 2^7, 2^8];
        for j = 1:length(NSteps)
            for RK = RKMethods
                nSteps = NSteps(j);
                dt = 1.0/nSteps;
                
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
                ThisInfo(RK,j).F = ThisInfo(RK,j).F(1:3:end);
                
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
            print(['ExampleThree-Prev-Residual-', ElementType], '-dpdf')
            figure(2106)
            pause(1)
            print(['ExampleThree-Prev-Bearing-', ElementType], '-dpdf')
            figure(2107)
            pause(1)
            print(['ExampleThree-Prev-Error-', ElementType], '-dpdf')
            figure(2108)
            pause(1)
            print(['ExampleThree-Prev-TimeError-', ElementType], '-dpdf')
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
