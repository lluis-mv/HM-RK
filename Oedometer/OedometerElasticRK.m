function [] = OedometerElasticRK()

addpath('../')
% 1. Define the problem

T = 1E-4;


CP.E = 1000;
CP.nu = 0.3;
CP.k = 1E-2;
nu = CP.nu;
CP.M = CP.E*(1-nu)/(1+nu)/(1-2*nu);

t = T/CP.M/CP.k;


eSize = 0.02;

model = createpde(1);

dx = 0.05; dy = 1;
R1 = [3,4,0, dx, dx, 0, 0, 0, dy, dy]';
g = decsg(R1);
geometryFromEdges(model, g);

mesh = generateMesh(model, 'Hmax', eSize, 'GeometricOrder','linear');

Nodes = mesh.Nodes';
Elements = mesh.Elements';

mesh = generateMesh(model, 'Hmax', eSize);

Nodes2 = mesh.Nodes';
Elements2 = mesh.Elements';




NSteps = 10.^linspace(0, 5, 8);
NSteps = 10.^linspace(0, 3.8, 8);
NSteps = 10.^[0:0.5:5];
NSteps = floor(NSteps); NSteps = sort(NSteps);



firstTime = true;

figure(99)
for i = 1:4; subplot(2,2,i); hold off; end


% for RKMethod = [1:8]
% for RKMethod = [1:8]
for RKMethod = 1:8    
    
    i = 1;
    
    for nSteps = NSteps
        
        dt = t/nSteps;
%         Nodes = Nodes2;
%         Elements = Elements2;
        [U,GPInfo] = ComputeThisLinearProblem(Nodes, Elements, CP, dt, nSteps, 'T3T3', RKMethod);
        
        % %         figure(904)
        % %         plot(U(3:3:end), Nodes(:,2), 'b*')
        % %
        % %         figure(905)
        % %         plot(U(2:3:end), Nodes(:,2), 'b*')
        
        
        if ( firstTime)
            [Xa] = ComputeAnalyticalSolution(Nodes, Elements, 'T3T3', t, CP, GPInfo,U);
            firstTime = false;
        end
        [L2(i), L2U(i), LInf(i), LInfU(i)] = ComputeErrorNorms(U, Xa, Nodes, Elements, GPInfo, CP);
        
        
        
        %         [U,GPInfo] = ComputeThisLinearProblem(Nodes2, Elements2, CP, dt, nSteps, 'T6T6');
        %         if ( firstTime)
        %             [Xa2] = ComputeAnalyticalSolution(Nodes2, Elements2,'T6T6',  t, CP, GPInfo,U);
        %         end
        %         [L2i(i), L2Ui(i), LInfi(i), LInfUi(i)] = ComputeErrorNorms(U, Xa2, Nodes2, Elements2, GPInfo, CP);
        %
        %         figure(904)
        %         hold on
        %         plot(U(3:3:end), Nodes2(:,2), 'r*')
        %
        %
        %         figure(905)
        %         hold on
        %         plot(U(2:3:end), Nodes2(:,2), 'r*')
        %
        %
        %
        %         [U,GPInfo] = ComputeThisLinearProblem(Nodes2, Elements2, CP, dt, nSteps, 'T6T3');
        %         if ( firstTime)
        %             [Xa3] = ComputeAnalyticalSolution(Nodes2, Elements2,'T6T3',  t, CP, GPInfo,U);
        %             firstTime = false;
        %         end
        %         [L2g(i), L2Ug(i), LInfg(i), LInfUg(i)] = ComputeErrorNorms(U, Xa3, Nodes2, Elements2, GPInfo, CP);
        %
        % %         figure(904)
        % %         hold on
        % %         plot(U(3:3:end), Nodes2(:,2), 'g*')
        % %         hold off
        % %
        % %         figure(905)
        % %         hold on
        % %         plot(U(2:3:end), Nodes2(:,2), 'g*')
        % %         hold off
        
        
        
        
        
        
        
        i = i+1;
        
    end
    figure(99)
    subplot(2,2,1)
    loglog( NSteps, L2, '*-.')
    xlabel('nSteps')
    ylabel('L2');
    hold on
    
    
    subplot(2,2,2)
    loglog( NSteps, LInf, '*-.')
    xlabel('nSteps')
    ylabel('LInf');
    hold on
    
    
    subplot(2,2,3)
    loglog( NSteps, L2U, '*-.')
    xlabel('nSteps')
    ylabel('L2 u');
    hold on
    
    
    subplot(2,2,4)
    loglog( NSteps, LInfU, '*-.')
    xlabel('nSteps')
    ylabel('LInf u');
    hold on
end




function [Xa] = ComputeAnalyticalSolution(Nodes, Elements, ElementType, t, CP, GPInfo, Xnum)
Xa = 0*Xnum;

% analytical solution
[Ca, Ka ] = EnsambleMatrices(Nodes, Elements, GPInfo, CP, ElementType, 3, 3, false, 0);

[Ca, Ka, X0, ~] = ApplyBoundaryConditions(Nodes, Elements, Ca, Ka);

Aa = Ca\(Ka);

[vectors, values] = eig(full(Aa), 'nobalance');

Xa = 0*Xnum;

c = (vectors)\X0;

for i = 1:size(values, 1)
    Xa = Xa + c(i)*exp(values(i,i)*t)*vectors(:,i);
end


Xa = real(Xa);



function [L2, L2U, LInf, LInfU] = ComputeErrorNorms(X, Xa, Nodes, Elements, GPInfo, CP)



nNodes = size(Nodes, 1);
nElements = size(Elements, 1);


indexWP = 3*[1:nNodes];
LInf = max( abs( X(indexWP)-Xa(indexWP)));
LInfU = 0;
for i = 1:nNodes
    ind = 3*(i-1)+[1,2];
    thisNorm = norm(Xa(ind)-X(ind));
    LInfU = max(LInfU, thisNorm);
end

L2 = 0;
L2U = 0;


alfa = 2/3; beta = 1/6;
N1 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 1/6;
N2 = [ 1 - alfa - beta, alfa,  beta];
alfa = 1/6; beta = 2/3;
N3 = [ 1 - alfa - beta, alfa,  beta];
% %
% % for el = 1:nElements
% %     Cel = Elements(el,:);
% %     indWP = 3*(Cel-1)+3;
% %     err = ( Xa(indWP)-X(indWP));
% %     L2 = L2 + GPInfo(el).Weight/3* ( abs(N1*err) + abs(N2*err)+abs(N3*err));
% %
% %     indx = 3*(Cel-1)+1;
% %     indy = 3*(Cel-1)+2;
% %     ux = Xa(indx)-X(indx);
% %     uy = Xa(indy)-X(indy);
% %     L2U = L2U + GPInfo(el).Weight/3* ( norm(N1*[ux,uy]) + norm(N2*[ux,uy])+norm(N3*[ux,uy]));
% %
% % end



LInfU = LInfU*CP.M;
L2U = L2U*CP.M;

