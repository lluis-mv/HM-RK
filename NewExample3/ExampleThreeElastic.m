function [] = ExampleThreeElastic()



addpath('../')
% 1. Define the problem


CP.HydroMechanical = true;
CP.E = 1000;
CP.nu = 0.3;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);
CP.k = 1E-8;
CP.Elastic = true;
CP.MCC = true;

ESIZE = 0.35;

RKReference = 8;
RKMethods = [8,1:7];



LinearElastic = false;
if ( LinearElastic)
    CP.Elastic = true;
    CP.MCC = false;
end

[NodesQ, ElementsQ] = ReadTheMesh('ThisMesh.msh');
[NodesT, ElementsT] = ConvertToTriangles(NodesQ, ElementsQ);


for Elem = [1:2]
    
    RKMethods = [8,1:7, 9];
    
    
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
        
        figure(1);
        PlotMesh(Nodes, Elements)
        drawnow
        axis equal
        axis off
        my_Print(['ExampleThree-Elastic-FemMesh-', ElementType], '-dpdf')
        
        % Estimate the element size
        
        
        Nadim = 20;
        
        NSteps = [4, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8];
        for j = 1:length(NSteps)
            for RK = RKMethods
                nSteps = NSteps(j);
                dt = 1.0/nSteps;
                
                disp(RK)
                if ( RK < 9)
                    tic;
                    
                    if ( LinearElastic)
                        [U, GPInfo,  information] = ComputeLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1);
                        rrr = 10E-14*rand();
                    else
                        [U,GPInfo, rrr, information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1, false);
                    end
                    timing = toc;
                elseif ( RK == 9)
                    tic;
                    [U, GPInfo, rrr,  information] = ComputeImplicitNonLinearProblem(Nodes, Elements, CP, dt, nSteps, ElementType);
                    timing = toc;
                    if ( isnan(rrr) )
                        timing = nan;
                        for mnm = 1:length(information)
                            for nmn = 1:length(information(mnm).F)
                                information(mnm).F(nmn)= nan;
                            end
                        end
                    end
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
                ind = find(N == 0);
                N(ind) = nan; ddtt(ind) = nan; Time(ind) = nan; RES(ind) = nan;
                
                if ( j < 2 )
                    continue;
                end
                
                figure(2105); clf
                for jj = 1:size(RES,1)
                    merda = '*-.';
                    if ( jj == 8)
                        merda = 'k*-.';
                    elseif (jj == 9)
                        merda = 'rs-.';
                    end
                    loglog(ddtt(jj,:), RES(jj,:), [merda])
                    hold on
                    
                end
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Residual', 'interpreter', 'latex')
                if (length(RKMethods) ==9)
                    ll = legend('1','2','3','4','5','6','7','8', 'I', 'location','bestoutside');
                else
                    ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                end
                set(ll, 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                figure(2106); clf
                for jj = 1:size(RES,1)
                    merda = '*-.';
                    if ( jj == 8)
                        merda = 'k*-.';
                    elseif (jj == 9)
                        merda = 'rs-.';
                    end
                    loglog(ddtt(jj,:), N(jj,:), [merda])
                    hold on
                end
                grid minor
                drawnow
                hold off
                if (length(RKMethods) ==9)
                    ll = legend('1','2','3','4','5','6','7','8', 'I', 'location','bestoutside');
                else
                    ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                end
                set(ll, 'interpreter', 'latex')
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Footing load', 'interpreter', 'latex')
                
                figure(2107); clf
                figure(2108); clf
                ind = find(N == 0);
                N(ind) = nan;
                ind = find(Time == 0);
                Time(ind) = nan;
                for jj = 1:size(RES,1)
                    merda = '*-.';
                    if ( jj == 8)
                        merda = 'k*-.';
                    elseif (jj == 9)
                        merda = 'rs-.';
                    end
                    plotThis = abs(N(jj,:)-Nadim);
                    index = find( plotThis == 0);
                    plotThis(index) = 1E-14*(1+rand(size(index)));
                    
                    figure(2107)
                    loglog(ddtt(jj,:), plotThis, [merda])
                    hold on
                    drawnow
                    
                    figure(2108)
                    loglog(Time(jj,:), plotThis, [merda])
                    hold on
                    drawnow
                    
                end
                figure(2107)
                if (length(RKMethods) ==9)
                    ll = legend('1','2','3','4','5','6','7','8', 'I', 'location','bestoutside');
                else
                    ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                end
                set(ll, 'interpreter', 'latex')
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                
                figure(2108)
                if (length(RKMethods) ==9)
                    ll = legend('1','2','3','4','5','6','7','8', 'I', 'location','bestoutside');
                else
                    ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                end
                set(ll, 'interpreter', 'latex')
                xlabel('Computational cost (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                
                % Printing....
                figure(2105)
                pause(1)
                drawnow
                my_Print(['ExampleThree-Elastic-Residual-', ElementType], '-dpdf')
                figure(2106)
                pause(1)
                drawnow
                my_Print(['ExampleThree-Elastic-Bearing-', ElementType], '-dpdf')
                figure(2107)
                pause(1)
                drawnow
                my_Print(['ExampleThree-Elastic-Error-', ElementType], '-dpdf')
                figure(2108)
                pause(1)
                drawnow
                my_Print(['ExampleThree-Elastic-TimeError-', ElementType], '-dpdf')
            end
        end
        clear Time;
        clear RES;
        clear ddtt;
        clear N;
        clear Time
    end
end


function [] = my_Print(XNAME, FORMAT)
print(XNAME, FORMAT)
savefig([XNAME, '.fig'])
