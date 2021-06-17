function [ProblemData] = SolveProblemAndComputeNorms( ProblemData)


if ( ProblemData.PrimalForm)
    if ( ProblemData.StabMethod < 7)
        ProblemData = ConstructAndSolveProblem( ProblemData );
    elseif ( ProblemData.StabMethod ==10)
        ProblemData = ConstructAndSolveProblemDiffOrder( ProblemData);
    elseif ( ProblemData.StabMethod == 11)
        ProblemData = ConstructAndSolveProblemDiffOrderP3P2( ProblemData);
    elseif ( ProblemData.StabMethod == 12)
        ProblemData = ConstructAndSolveProblemDiffOrderP4P3( ProblemData);
    elseif ( ProblemData.StabMethod ==14)
        ProblemData = ConstructAndSolveProblemDiffOrderP2P2( ProblemData);
    elseif ( ProblemData.StabMethod == 15)
        ProblemData = ConstructAndSolveProblemDiffOrderP3P1( ProblemData);
    elseif ( ProblemData.StabMethod == 16)
        ProblemData = ConstructAndSolveProblemDiffOrderP4P1( ProblemData);
    
    else
        error 'no problem defined';
    end
else
    ProblemData = ConstructAndSolveProblemUWwP( ProblemData );
end

TT = ProblemData.TT;


% L2 water pressure
xxGP = GetGaussPointPosition(ProblemData.XX);
WaterPressureGPh = EvaluateNodalAtGP(ProblemData.NodalWaterPressure, ProblemData.OrderPW);
WaterPressureGPa = EvalConsolidation(xxGP, TT);
L2 = EvaluateL2Norm( ProblemData.XX(2)-ProblemData.XX(1), WaterPressureGPh, WaterPressureGPa);

% LInf water pressure 
xx = ProblemData.XX;
if ( isempty(ProblemData.xxPW) == false)
    xx = ProblemData.xxPW;
end
wpNodes = EvalConsolidation(xx, TT);
LInf = max( abs(ProblemData.NodalWaterPressure - wpNodes));

% L2Displ
dispGPh = EvaluateNodalAtGP(ProblemData.NodalDisplacement, ProblemData.OrderU);
dispGPa = EvalConsolidationDISPL(xxGP, TT)/ProblemData.M;
L2DISPL = EvaluateL2Norm( ProblemData.XX(2)-ProblemData.XX(1), dispGPh, dispGPa);

% LInf displ
xx = ProblemData.XX;
if ( isempty(ProblemData.xxU) == false)
    xx = ProblemData.xxU;
end
uNodes = EvalConsolidationDISPL(xx, TT)/ProblemData.M;
LInfDISPL = max( abs(ProblemData.NodalDisplacement - uNodes));

if ( ProblemData.PrimalForm == false)
    dispWGPh = EvaluateNodalAtGP(ProblemData.NodalWaterDisplacement, ProblemData.OrderU);
    L2WATER = EvaluateL2Norm( ProblemData.XX(2)-ProblemData.XX(1), dispWGPh, -dispGPa);
    
    % LInf displ
    xx = ProblemData.XX;
    if ( isempty(ProblemData.xxU) == false)
        xx = ProblemData.xxU;
    end
    uWNodes = -EvalConsolidationDISPL(xx, TT)/ProblemData.M;
    LInfWATER = max( abs(ProblemData.NodalWaterDisplacement - uWNodes));

else
    LInfWATER = 0;
    L2WATER = 0;
end

Proportion = ProblemData.M;
Proportion = 1;

ProblemData.LInfDISP = LInfDISPL*Proportion;
ProblemData.L2DISP = L2DISPL*Proportion;

ProblemData.LInfWATER = LInfWATER*Proportion;
ProblemData.L2WATER = L2WATER*Proportion;

ProblemData.LInf = LInf;
ProblemData.L2 = L2;


figure(101)
if ( isempty(ProblemData.xxU))
    xx = ProblemData.XX;
    plot(xx, ProblemData.NodalDisplacement, '*-', xx, uNodes)
else
    plot(xx, ProblemData.NodalDisplacement, '*k', xx, uNodes)
    hold on
    [xx,yy] = EvaluateSolutionBeatiful(ProblemData.XX, ProblemData.NodalDisplacement, ProblemData.OrderU);
    plot(xx, yy, 'k')
    hold off
end
    

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

if ( isempty(ProblemData.NodalWaterDisplacement) == false)
    if ( length(xx) == length(ProblemData.NodalWaterDisplacement) )
        figure(105)
        plot(xx, ProblemData.NodalWaterDisplacement, '*-', xx, -uNodes)
    end
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
        FieldGP(i,:) = N'*Field(i:i+1);
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




function pw = EvalConsolidation(XX, TT)
pw = 0*XX;
nTermsZero = 0;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:5000
            aux = pi/2*(2*m+1);
            term = 2/aux * sin( aux * XX(i,j)) * exp( - aux^2 * TT);
            if ( abs(term) < 1e-12)
                nTermsZero = nTermsZero + 1;
            else
                nTermsZero = 0;
            end
            pw(i,j) = pw(i,j)+term ;
            if ( nTermsZero > 20)
                break;
            end
        end
    end
end

function uu = EvalConsolidationDISPL(XX, TT)
uu = -(1-XX);
nTermsZero = 0;
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:5000
            z=XX(i,j);
            term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
            if ( abs(term) < 1e-12)
                nTermsZero = nTermsZero + 1;
            else
                nTermsZero = 0;
            end
            uu(i,j) = uu(i,j)+term;
            if ( nTermsZero > 20)
                break;
            end
        end
    end
end



