
function uu = EvalConsolidationDISPL(XX, TT)
uu = -(1-XX);
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        for m = 0:100
            z=XX(i,j);
            term = +(exp(-(TT*pi^2*(2*m + 1)^2)/4)*(8*sin(pi*m) + 8*cos((z*pi*(2*m + 1))/2)))/(pi^2*(2*m + 1)^2);
            uu(i,j) = uu(i,j)+term;
        end
    end
end
