function [] = ExampleThree()



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
CP.RK = 1;

ESIZE = 0.35;

RKReference = 8;
RKMethods = [8,1:7];

RKReference = 8;
for Elem = [1,2,3]
    RKMethods = [8,1:7];
    if ( Elem == 2)
        RKMethods = [8,1:7, 9];
    end
    
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
        print('ExampleThree-Plastic-FemMesh', '-dpdf')
        
        % Estimate the element size
        
        mesha = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');
        Nodesa = mesha.Nodes';
        Elementsa = mesha.Elements';
        [GPInfo] = ComputeElementalMatrices(Nodesa, Elementsa, CP, 'T3T3');
        he = mean(sqrt( mean([GPInfo(:,:).Weight])));
        
        Nadim = 20;
        
        NSteps = [2^5, 2^6, 2^7, 2^8, 2^9];
        for j = 1:length(NSteps)
            for RK = RKMethods
                nSteps = NSteps(j);
                dt = 0.15/nSteps;
                
                disp(RK)
                if ( RK < 9)
                	tic;
	                [U,GPInfo, rrr, information] = ComputeNLProblem(Nodes, Elements, CP, dt, nSteps, ElementType, RK, 1);
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
                ThisInfo(RK,j).F = ThisInfo(RK,j).F(1:2:end);
                
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
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    elseif (jj == 9)
                        merda = 'r';
                    end
                    loglog(ddtt(jj,:), RES(jj,:), [merda, '*-.'])
                    hold on
                    
                end
                xticks([4E-4, 1E-3, 4E-3])
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
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    elseif (jj == 9)
                        merda = 'r';
                    end
                    loglog(ddtt(jj,:), N(jj,:), [merda, '*-.'])
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
                xticks([4E-4, 1E-3, 4E-3])
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Footing load', 'interpreter', 'latex')
                
                figure(2107); clf
                ind = find(N == 0);
                N(ind) = nan;
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    elseif (jj == 9)
                        merda = 'r';
                    end
                    plotThis = abs(N(jj,:)-Nadim);
                    if ( jj ~= RKReference)
                        index = find( plotThis == 0);
                        plotThis(index) = 1E-14*(1+rand(size(index)));
                    end
                    loglog(ddtt(jj,:), plotThis, [merda, '*-.'])
                    hold on
                end
                
                if (length(RKMethods) ==9)
                    ll = legend('1','2','3','4','5','6','7','8', 'I', 'location','bestoutside');
                else
                ll = legend('1','2','3','4','5','6','7','8', 'location','bestoutside');
                end
                set(ll, 'interpreter', 'latex')
                xticks([4E-4, 1E-3, 4E-3])
                xlabel('$\Delta t$ (s)', 'interpreter', 'latex')
                ylabel('Error norm', 'interpreter', 'latex')
                grid minor
                drawnow
                hold off
                
                
                figure(2108); clf
                ind = find(Time == 0);
                Time(ind) = nan;
                for jj = 1:size(RES,1)
                    merda = '';
                    if ( jj == 8)
                        merda = 'k';
                    elseif (jj == 9)
                        merda = 'r';
                    end
                    plotThis = abs(N(jj,:)-Nadim);
                    if ( jj ~= RKReference)
                        index = find( plotThis == 0);
                        plotThis(index) = 1E-14*(1+rand(size(index)));
                    end
                    loglog(Time(jj,:), plotThis, [merda, '*-.'])
                    hold on
                end
                
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
                print(['ExampleThree-Plastic-Residual-', num2str(Elem)], '-dpdf')
                figure(2106)
                pause(1)
                print(['ExampleThree-Plastic-Bearing-', num2str(Elem)], '-dpdf')
                figure(2107)
                pause(1)
                print(['ExampleThree-Plastic-Error-', num2str(Elem)], '-dpdf')
                figure(2108)
                pause(1)
                print(['ExampleThree-Plastic-TimeError-', num2str(Elem)], '-dpdf')
            end
            
            
            
        end
        clear Time;
        clear RES;
        clear ddtt;
        clear N;
        clear Time
    end
end



