

x0 = [ 2/5, -1/5, -1/5, 3/5, -1/5, 3/5 -3/5, 3/5, 0, 0, 4/5, -4/5 -3/5, 0, 3/5, -4/5, 4/5, 0]';
a = x0(1:6);
b = x0(7:12);
c = x0(13:18);

p = [-1,0,0,0,0.0,0.0];


ALFA = linspace(0,1,50);
BETA = linspace(0,1,50);
[ALFA, BETA] = meshgrid(ALFA, BETA);
T = delaunay(ALFA,BETA);
%trimesh(T,x,y,z)
for i = 1:size(ALFA,1)
    for j = 1:size(ALFA,2)
        alfa = ALFA(i,j);
        beta = BETA(i,j);
        N =  [ (1 - alfa - beta)*(1-2*alfa-2*beta);
            alfa*(2*alfa-1);
            beta*(2*beta-1);
            4*alfa*(1-alfa-beta);
            4*alfa*beta;
            4*beta*(1-alfa-beta)];
        
        na = a+b*alfa+c*beta;
        ph(i,j) = p*N;
        p2(i,j) = p*na;
 
    end
end

[rs] = find( BETA > 1-ALFA -0.001);
ph(rs) = nan;
p2(rs) =nan;
figure(232)
trimesh(T, ALFA, BETA, ph);
hold on
trimesh(T, ALFA, BETA, p2);
hold off