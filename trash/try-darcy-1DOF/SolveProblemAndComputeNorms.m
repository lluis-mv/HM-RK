function [ProblemData] = SolveProblemAndComputeNorms( ProblemData)


if ( ProblemData.PrimalForm)
    ProblemData = ConstructAndSolveProblem( ProblemData );    
end
TT = ProblemData.tFinal;


nGP = 5;

% L2Displ
xxGP = GetGaussPointPosition(ProblemData.XX);
dispGPh = EvaluateNodalAtGP(ProblemData.NodalDisplacement, ProblemData.OrderU);
dispGPa = EvalSolution(xxGP, TT);
L2DISPL = EvaluateL2Norm( ProblemData.XX(2)-ProblemData.XX(1), dispGPh, dispGPa);

% LInf displ
xx = ProblemData.XX;
if ( isempty(ProblemData.xxU) == false)
    xx = ProblemData.xxU;
end
uNodes = EvalSolution(xx, TT);
LInfDISPL = max( abs(ProblemData.NodalDisplacement - uNodes));

ProblemData.LInfDISP = LInfDISPL;
ProblemData.L2DISP = L2DISPL;




figure(101)
if ( isempty(ProblemData.xxU) || ProblemData.OrderU == 1)
    xx = ProblemData.XX;
    plot(xx, ProblemData.NodalDisplacement, '*-', xx, uNodes)
else
    plot(xx, ProblemData.NodalDisplacement, '*k', xx, uNodes)
    hold on
    [xx,yy] = EvaluateSolutionBeatiful(ProblemData.XX, ProblemData.NodalDisplacement, ProblemData.OrderU);
    plot(xx, yy, 'k')
    hold off
end
pause(0.01)    
return;
figure(102)
if ( isempty(ProblemData.xxPW) )
    xx = ProblemData.XX;
    plot(xx, ProblemData.NodalWaterPressure, '*-', xx, wpNodes)
else
    xx = ProblemData.xxPW;
    plot(xx, ProblemData.NodalWaterPressure, '*k', xx, wpNodes)
    hold on
    [xx1,yy1] = EvaluateSolutionBeatiful(ProblemData.XX, ProblemData.NodalWaterPressure, ProblemData.OrderPW);
    plot(xx1, yy1, 'k')
    hold off
end

pause(0.0001)


return
figure(103)
xPlot = [];
yPlot = [];
zPlot = [];
for a = 1:size(xxGP, 1)
    xPlot = [xPlot, xxGP(a,:)];
    yPlot = [yPlot, WaterPressureGPh(a,:)];
    zPlot = [zPlot, WaterPressureGPa(a,:)];
end
plot(xPlot, yPlot, 'r', xPlot, zPlot, 'g')
hold on
plot(xx, ProblemData.NodalWaterPressure, '*k', xx, wpNodes)
legend('GP-num', 'GP-an', 'Nodal-h', 'Nodal-an', 'location','best')
hold off

xx = ProblemData.xxU;
if ( isempty(ProblemData.xxU))
    xx = ProblemData.XX;
end
figure(104)
xPlot = [];
yPlot = [];
zPlot = [];
for a = 1:size(xxGP, 1)
    xPlot = [xPlot, xxGP(a,:)];
    yPlot = [yPlot, dispGPh(a,:)];
    zPlot = [zPlot, dispGPa(a,:)];
end
plot(xPlot, yPlot, 'r', xPlot, zPlot, 'g')
hold on
plot(xx, ProblemData.NodalDisplacement, '*k', xx, uNodes)
legend('GP-num', 'GP-an', 'Nodal-h', 'Nodal-an', 'location','best')
hold off


function xxGP = GetGaussPointPosition(XX)

chis = [ -1/3*sqrt(5+2*sqrt(10/7)), -1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];

SizeField = length(XX);

N = [-(chis-1)/2; (chis+1)/2];
xxGP = zeros(SizeField-1,5);
for i = 1:SizeField-1
    xxGP(i,:) = N'*XX(i:i+1);
end
        
        


function FieldGP = EvaluateNodalAtGP(Field, Order)

chi = [ -1/3*sqrt(5+2*sqrt(10/7)), -1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];

SizeField = length(Field);

if ( Order == 1)
    N = [-(chi-1)/2; (chi+1)/2];
    FieldGP = zeros(SizeField-1,5);
    for i = 1:SizeField-1
        FieldGP(i,:) = N'*Field(i:i+1)';
    end
elseif ( Order == 2)
    SizeField = (SizeField-1)/2;
    N = [1/2*chi.*(chi-1); (1-chi.^2); 1/2*chi.*(chi+1)];
    FieldGP = zeros(SizeField-1,5);
    j = 1;
    for i = 1:1:SizeField
        FieldGP(i,:) = N'*(Field(j:j+2)');
        j = j+2;
    end
elseif( Order == 3)
    N = [ -9/16*(chi+1/3).*(chi-1/3).*(chi-1); %N1
    27/16*(chi+1).*(chi-1/3).*(chi-1); %N2
    -27/16*(chi+1).*(chi+1/3).*(chi-1); %N3
    9/16*(chi+1).*(chi+1/3).*(chi-1/3) %N4
    ];
    SizeField = (SizeField-1)/3;
    FieldGP = zeros(SizeField-1,5);
    j = 1;
    for i = 1:1:SizeField
        FieldGP(i,:) = N'*(Field(j:j+3)');
        j = j+3;
    end
elseif ( Order == 4)
    N = [    (2*chi.*(chi - 1).*(chi - 1/2).*(chi + 1/2))/3;
    -(8*chi.*(chi - 1).*(chi + 1).*(chi - 1/2))/3;
 4*(chi - 1).*(chi + 1).*(chi - 1/2).*(chi + 1/2);
    -(8*chi.*(chi - 1).*(chi + 1).*(chi + 1/2))/3;
   (2*chi.*(chi + 1).*(chi - 1/2).*(chi + 1/2))/3;
   ];
    SizeField = (SizeField-1)/4;
    FieldGP = zeros(SizeField-1,5);
    j = 1;
    for i = 1:1:SizeField
        FieldGP(i,:) = N'*(Field(j:j+4)');
        j = j+4;
    end
else
    error;
end
        
function [xx,yy] = EvaluateSolutionBeatiful(XX, UU, Order)        
xx=[]; yy=[];
chi=linspace(-1,1);
SizeField = length(UU);
Nx = [-(chi-1)/2; (chi+1)/2];
if (Order == 2)
    N = [1/2*chi.*(chi-1); (1-chi.^2); 1/2*chi.*(chi+1)];
    toAdd = 2;
    SizeField = (SizeField-1)/2;
elseif (Order == 3)
     N = [ -9/16*(chi+1/3).*(chi-1/3).*(chi-1); %N1
    27/16*(chi+1).*(chi-1/3).*(chi-1); %N2
    -27/16*(chi+1).*(chi+1/3).*(chi-1); %N3
    9/16*(chi+1).*(chi+1/3).*(chi-1/3) %N4
    ];
    toAdd = 3;
    SizeField = (SizeField-1)/3;
elseif ( Order == 4)
      N = [    (2*chi.*(chi - 1).*(chi - 1/2).*(chi + 1/2))/3;
    -(8*chi.*(chi - 1).*(chi + 1).*(chi - 1/2))/3;
 4*(chi - 1).*(chi + 1).*(chi - 1/2).*(chi + 1/2);
    -(8*chi.*(chi - 1).*(chi + 1).*(chi + 1/2))/3;
   (2*chi.*(chi + 1).*(chi - 1/2).*(chi + 1/2))/3;
   ];
toAdd = 4;
    SizeField = (SizeField-1)/4;
end
j=1;
for i = 1:1:SizeField
    xx = [xx; Nx'*XX(i:i+1)];
    yy = [yy; N'*(UU(j:j+toAdd)')];
    j = j+toAdd;
end

function L2 = EvaluateL2Norm( h, Num, Anal)

L2 = 0;

w = (h/2)*[(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225, (322+13*sqrt(70))/900,(322-13*sqrt(70))/900]';

for i = 1:size(Num,1)
    L2 = L2 + (   (Num(i,:)-Anal(i,:)).^2 * w);
end
L2 = sqrt(L2);




function uu = EvalSolution(XX, TT)
uu = solution(XX, TT);
return;
